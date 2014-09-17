#!/usr/bin/python
# this program is used to classify a certain number of reads that are
# in top K best paths.
import os
import sys
import operator
import subprocess
import networkx as nx
from pprint import pprint
import random
import heapq
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein
            

# define the read class.
class Read:
    def __init__(self, name):
        self.name = name
        self.seq = ''
        self.members = set()
  
    def __repr__(self):
        return repr((self.name, self.score, 
                     self.begin_state, self.end_state, 
                     self.seq, self.members))
   
    def set_seq(self, seq):
        self.seq = seq
        self.length = len(seq)
   
    def set_begin_state(self, begin_state):
        self.begin_state = begin_state

    def set_end_state(self, end_state):
        self.end_state = end_state
 
    def set_score(self, score):
        self.score = score

    def add_member(self, member):
        self.members.add(member)

# get domains with reads aligned to them.
def get_aligned_read_dict(in_file_name, target_domain):
    aligned_read_dict = {} # reads that belong to domains.
    with open(in_file_name, 'Ur') as f:
        for line in f:
            if not line.strip():
                continue
            items = line.strip().split()
            domain = items[1].split('=')[-1][:7]
            if domain == target_domain:
                read_name = items[0][1:]
                score = float(items[2].split('=')[-1])
                begin_state = int(items[4].split('=')[-1])
                end_state = int(items[5].split('=')[-1])
                strand_name = items[-1].split('=')[1]
                if strand_name == 'plus':
                    strand = '+'
                else:
                    strand = '-'
                composite_read_name = read_name + '$' + strand 
                if composite_read_name in aligned_read_dict and \
                    score <= aligned_read_dict[composite_read_name].score:
                    continue
                read = Read(composite_read_name)
                read.set_score(score)
                read.set_begin_state(begin_state)
                read.set_end_state(end_state)
                read.add_member(read_name)
                aligned_read_dict[composite_read_name] = read
    return aligned_read_dict

# remove reads that are aligned to the same state with lower scores.
def get_trimmed_aligned_read_dict(aligned_read_dict):
    begin_state_dict = {}
    for read in aligned_read_dict.values(): 
        begin_state = read.begin_state
        score = read.score
        begin_state_dict.setdefault(begin_state, read)
        if begin_state_dict[begin_state].score < score:
            begin_state_dict[begin_state] = read
    trimmed_aligned_read_dict = {}
    for read in begin_state_dict.values():
        trimmed_aligned_read_dict[read.name] = read
    return trimmed_aligned_read_dict 

# get the sequence of the reads.
def set_read_seq(fasta_file_name, read_dict):
    with open(fasta_file_name, 'Ur') as f:
        for record in SeqIO.parse(f, 'fasta'):
            read_name = record.id
            composite_read_name = read_name + '$+'
            if composite_read_name in read_dict:
                seq = str(record.seq)
                read_dict[composite_read_name].set_seq(seq)
            composite_read_name = read_name + '$-'
            if composite_read_name in read_dict:
                seq = str(record.seq.reverse_complement())
                read_dict[composite_read_name].set_seq(seq)

def get_compressed_read_dict(target_read_dict, prefix):
    seq_read_dict = {}
    for read_name in target_read_dict:
        seq = target_read_dict[read_name].seq
        seq_read_dict.setdefault(seq, [])
        seq_read_dict[seq].append(read_name)
    compressed_read_dict = {}
    index = 1
    for seq in seq_read_dict:
        tag_name = '%s%d' % (prefix, index)
        index += 1
        # use the first read that has the same sequence.
        first_read_name = seq_read_dict[seq][0]
        first_read = target_read_dict[first_read_name]
        tag = Read(tag_name)
        tag.set_begin_state(first_read.begin_state)
        tag.set_end_state(first_read.end_state)
        tag.set_score(first_read.score)
        tag.set_seq(first_read.seq)
        # keep the information
        for read_name in seq_read_dict[seq]:
            tag.add_member(read_name)
        compressed_read_dict[tag_name] = tag
    return compressed_read_dict

# get the set of reads mapped to the target domain.
def get_mapped_read_set(in_file_name, target_domain):
    mapped_read_set = set() 
    with open(in_file_name, 'Ur') as f:
        for line in f:
            if not line.strip():
                continue
            items = line.rstrip().split()
            domain = items[0][0:7]
            if domain == target_domain:
                read_name = items[1]
                mapped_read_set.add(read_name)
    return mapped_read_set

# get the overlap length of two reads. return 0 if no overlap.
def get_seq_overlap_length(seq1, seq2):
    # no string should contain the other.
    hamming_thres = 4
    len1 = len(seq1)
    len2 = len(seq2)
    max_overlap = 0 # currently maximum overlap.
    # loop over number of overlappd pos.
    for i in xrange(min(len1, len2), 0, -1):
        if get_hamming_distance(seq1[-i:], seq2[:i]) <= hamming_thres:
            max_overlap = i
            break
    return max_overlap   

def get_hamming_distance(seq1, seq2):
    assert len(seq1) == len(seq2)
    hamming_distance = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            hamming_distance += 1
    return hamming_distance

# get the overlap of two positions.
# if there is no overlap, return negative value.
def get_pos_overlap_length(begin1, end1, begin2, end2):
    return min(end1, end2) - max(begin1, begin2) + 1

# add root to transfer node score to edge weight.
def add_root_to_subgraph(subgraph):
    for node in subgraph.nodes():
        if subgraph.in_degree(node) == 0:
            subgraph.add_edge('root', node, 
                              weight=subgraph.node[node]['score'],
                              overlap=0)
    subgraph.node['root']['type'] = 'negative'

# convert weights of reads into read number.
def get_read_num_from_set(read_set):
    read_num = 0
    for read in read_set:
        # root is not counted as a read.
        if read == 'root':
            continue
        read_weight = int(read.split('_')[-1])
        read_num += read_weight
    return read_num

# get confusion matrix.
def get_confusion_mat(mapped_reads, predicted_read_set, TEST_READ_NUM):
    mapped_num = get_read_num_from_set(mapped_reads)
    TP = get_read_num_from_set(mapped_reads & predicted_read_set)  
    FN = get_read_num_from_set(mapped_reads - predicted_read_set)
    FP = get_read_num_from_set(predicted_read_set - mapped_reads)
    TN = TEST_READ_NUM - mapped_num - FP
    return (TP, FN, FP, TN)


# get labels of nodes to display in graph.
# the current label is score of node.
def get_node_labels(subgraph):
    node_labels = {}
    node_labels['root'] = 'root'
    node_labels['positive'] = {}
    node_labels['negative'] = {}
    for read_name in subgraph.nodes():
        if read_name != 'root':
            score = subgraph.node[read_name]['score']
            if subgraph.node[read_name]['type']:
                node_labels['positive'][read_name] = '%.2f' % score
            else:
                node_labels['negative'][read_name] = '%.2f' % score
    return node_labels 

# get labels for edges.
# the current label is overlap between two nodes of the edge.
def get_edge_labels(subgraph):
    edge_labels = {}
    for node1, node2 in subgraph.edges():
        overlap = subgraph[node1][node2]['overlap']
        edge_labels[(node1, node2)] = overlap
    return edge_labels

# draw one graph. split subgraphs into different files.
def draw_graphs(G, folder_name):
    domain_name = G.graph['domain']
    dir = folder_name + '/' + domain_name
    if not os.path.exists(dir):
        os.makedirs(dir)
    subgraphs = nx.weakly_connected_component_subgraphs(G)
    add_root_to_subgraphs(subgraphs)
    subgraphs.sort(key=lambda subgraph: subgraph.number_of_nodes())
    for i in xrange(len(subgraphs)): 
        subgraph = subgraphs[i]
        pos = nx.spring_layout(subgraph)
        node_labels = get_node_labels(subgraph)
        positive_nodes = node_labels['positive'].keys()
        negative_nodes = node_labels['negative'].keys()
        labels = dict(node_labels['positive'], **(node_labels['negative']))
        edge_labels = get_edge_labels(subgraph)
        pl.figure(figsize=(16, 12))
        nx.draw_networkx_nodes(subgraph, pos, positive_nodes, alpha=0.5, node_color='w')
        nx.draw_networkx_nodes(subgraph, pos, negative_nodes, alpha=0.5, node_color='b')
        nx.draw_networkx_nodes(subgraph, pos, ['root'], node_color='g')
        nx.draw_networkx_edges(subgraph, pos, color='k')
        nx.draw_networkx_labels(subgraph, pos, labels, font_size=20)
        nx.draw_networkx_edge_labels(subgraph, pos, edge_labels, font_size=20)
        pl.axis('off')
        pl.savefig('%s/%s_subgraph_%d.png' % (dir, domain_name, i+1))

# draw the graph of a path. use red for positive reads 
# white for border reads and blue for negative reads.
def draw_path_subgraphs(paths, G):
    domain = G.graph['domain']
    folder_name = 'Graph/' + domain 
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    TP_list = []
    for i in range(len(paths)):
        path = paths[i]
        path_graph = G.subgraph(path)
        pos = nx.spring_layout(path_graph)
        contig = get_contig_from_path(path, G)
        contig_length = len(contig)
        file_name = '%s/contig%d_length%d.png' % (folder_name, i+1, contig_length)
        pl.figure(figsize=(16, 12))
        nodes = path_graph.nodes()
        node_colors = get_node_colors(nodes, path_graph)
        TP_list += filter(lambda x: G.node[x]['type']=='positive', nodes)
        edge_labels = get_edge_labels(path_graph)
        nx.draw_networkx_nodes(path_graph, pos, nodes, alpha=0.5, node_color=node_colors) 
        nx.draw_networkx_edges(path_graph, pos, color='k')
        nx.draw_networkx_labels(path_graph, pos, font_size=20)
        nx.draw_networkx_edge_labels(path_graph, pos, edge_labels, font_size=20)
        pl.axis('off')
        pl.savefig(file_name)


# return a list of colors corresponding to the same position of nodes.
def get_node_colors(nodes, subgraph):
    node_colors = []
    for node in nodes:
        if subgraph.node[node]['type'] == 'positive':
            node_colors.append('r') 
        elif subgraph.node[node]['type'] == 'border':
            node_colors.append('w')
        else: 
            node_colors.append('b')
    return node_colors

# output positive and negative node numbers of each subgraph for a domain.
def get_subgraph_size_list(subgraphs):
    subgraph_size_list = []
    for subgraph in subgraphs:
        num_pair = get_positive_negative_node_num(subgraph)
        subgraph_size_list.append(num_pair)
    subgraph_size_list.sort(key=lambda num_pair: sum(num_pair))
    return subgraph_size_list

# output information of the graph.
def output_graph_stat(G, mapped_read_lookup_dict):  
    domain = G.graph['domain']                  
    subgraphs = nx.weakly_connected_component_subgraphs(G)
    subgraph_num = len(subgraphs)
    subgraph_size_list = get_subgraph_size_list(subgraphs)
    mapped_read_num = len(mapped_read_lookup_dict[domain].keys())
    aligned_read_num = G.number_of_nodes()
    sys.stdout.write('%s:%d:%d:%d' % (domain, mapped_read_num, 
                                      aligned_read_num, subgraph_num))
    for positive_num, negative_num in subgraph_size_list:
        sys.stdout.write(' %d:%d' % (positive_num, negative_num)) 
    sys.stdout.write('\n')

# output the actual reads including meta information for subgraph.
def output_subgraph_data(subgraph, index, domain_read_dict, f):
    domain = subgraph.graph['domain']
    read_list = subgraph.nodes()
    read_list.sort(key=lambda read_name: subgraph.node[read_name]['begin_state'])
    for read_name in read_list:
        if subgraph.node[read_name]['type']:
            symbol = '+'
        else:
            symbol = '-'
        score = subgraph.node[read_name]['score']
        begin_state = subgraph.node[read_name]['begin_state']
        in_degree = subgraph.in_degree(read_name)
        out_degree = subgraph.out_degree(read_name)
        seq = domain_read_dict[domain][read_name].seq
        f.write('%s %s %d %s %.2f %d %d %d %s\n' % 
                (domain, read_name, index, symbol, score, begin_state, 
                 in_degree, out_degree, seq))

# output the actual reads including meta information for the domain.
def output_domain_data(G, domain_read_dict, folder_name):
    domain = G.graph['domain']
    dir = folder_name + '/' + domain
    if not os.path.exists(dir):
        os.makedirs(dir)
    out_file_name = dir + '/' + domain + '.data'
    subgraphs = nx.weakly_connected_component_subgraphs(G)
    subgraphs.sort(key=lambda subgraph:subgraph.number_of_nodes())
    with open(out_file_name, 'w') as f:
        for i in xrange(len(subgraphs)):
            subgraph = subgraphs[i]
            output_subgraph_data(subgraph, i+1, domain_read_dict, f)

# get the weighted length of a path.
def get_path_weighted_length(subgraph, path):
    assert path[0] == 'root'
    weighted_length = 0.0
    if path:
        for i in xrange(len(path)-1):
            weighted_length += subgraph[path[i]][path[i+1]]['weight']
    return weighted_length

# get the contig from a path.
def get_contig_from_path(path, G):
    contig = ''
    for i in xrange(1, len(path)):
        begin_node = path[i-1]
        end_node = path[i]
        overlap = G[begin_node][end_node]['overlap']
        seq = G.node[end_node]['seq']
        contig += seq[overlap:]
    return contig
            
# return top contigs and output them.
def output_top_N_contigs(k_longest_paths, G, out_file_name):
    with open(out_file_name, 'wb') as f:
        for i in range(len(k_longest_paths)):
            path = k_longest_paths[i]
            contig = get_contig_from_path(path, G)
            contig_length = len(contig)
            f.write('>contig%d_%d$%s' % 
                (i+1, contig_length, ':'.join(path[1:])+'\n')) # root is not written.
            f.write(contig+'\n') 

# classify N longest paths using hmmer.
# here contig means assembled long contigs that have hmmer hits.
def classify_top_N_contigs_with_hmmer(G, N, EVALUE_THRES, working_dir):
    if G.number_of_nodes() == 0:
        return set()
    K = N
    domain = G.graph['domain']
    k_longest_paths = get_k_longest_paths(G, K)
    dir = working_dir
    if not os.path.exists(dir):
        os.makedirs(dir)
    domain_hmmer_file_name = '%s/%s_top%d.hmmer' % (dir, domain, K)
    if os.path.exists(domain_hmmer_file_name):
        os.remove(domain_hmmer_file_name)
    generate_hmmer_output_file(G, k_longest_paths, dir, domain_hmmer_file_name, K, working_dir)
    classified_tag_set =  classify_hmmer_output_file(domain_hmmer_file_name,
                                                     EVALUE_THRES)
    return classified_tag_set

def get_classified_read_set(classified_tag_set, compressed_read_dict):
    classified_read_set = set()
    for tag in classified_tag_set:
        classified_read_set |= compressed_read_dict[tag].members
    return classified_read_set

def generate_hmmer_output_file(G, k_longest_paths, dir, domain_hmmer_file_name, K, working_dir):
    domain = G.graph['domain']
    #draw_path_subgraphs(k_longest_paths, G)
    current_path = os.path.dirname(os.path.realpath(__file__))
    hmmer_pipeline_file_name = '%s/hmmer3_pipeline.sh' % (current_path,)
    hmm_file_name = working_dir.rstrip('/') + '/HMMs/' + domain + '.hmm'
    for i in range(len(k_longest_paths)):
        path = k_longest_paths[i]
        contig_seq = Seq(get_contig_from_path(path, G), generic_dna)     
        contig_name = 'contig%d_%d$%s' % (i+1, len(contig_seq), ':'.join(path))
        contig_file_name = '%s/%s_top%d_contig%d.fna' % \
                            (dir, domain, K, i+1)
        #contig_seq_record = SeqRecord(contig_seq, id=contig_name, description='')
        #SeqIO.write(contig_seq_record, contig_file_name, 'fasta')
        for j in range(3):        
            frame_contig_seq = contig_seq[j:]
            peptide_seq = frame_contig_seq.translate(to_stop=False, stop_symbol="X")
            if len(peptide_seq) == 0:
                continue
            # in frame stop codon.
            #if len(peptide_seq) < len(contig_seq) / 3 - 1:
            #    continue
            peptide_file_name = '%s/%s_top%d_contig%d.frame%d' % \
                                 (dir, domain, K, i+1, j+1)
            peptide_seq_record = SeqRecord(peptide_seq, id=contig_name, description='')
            SeqIO.write(peptide_seq_record, peptide_file_name, 'fasta')
            peptide_hmmer_file_name = peptide_file_name + '.hmmer'
            # generate hmmer results.
            hmmer_cmd = '%s %s %s >%s' % (hmmer_pipeline_file_name, 
                                          hmm_file_name,
                                          peptide_file_name,
                                          peptide_hmmer_file_name)        
            subprocess.call(hmmer_cmd, shell=True)
            os.remove(peptide_file_name)
            # append the valid hmmer results into domain hmmer file.
            append_hmmer_cmd = 'cat %s >> %s' % (peptide_hmmer_file_name, domain_hmmer_file_name)
            subprocess.call(append_hmmer_cmd, shell=True)
            #os.remove(peptide_hmmer_file_name)

# classify reads in hmmer files.
def classify_hmmer_output_file(in_file_name, evalue_thres):
    classified_read_set = set()
    if not os.path.exists(in_file_name):
        return classified_read_set
    with open(in_file_name, 'Ur') as f:
        for line in f:
            if not line.strip():
                continue
            items = line.rstrip().split()
            contig_name = items[0]
            contig_length = int(items[0].split('$')[0].split('_')[1])
            read_list = items[0].split('$')[1].split(':')
            score = float(items[2])
            evalue = float(items[3])
            seq_begin = (int(items[6]) - 1) * 3 + 1
            seq_end = int(items[7]) * 3
            alignment_seq_length = seq_end - seq_begin + 1
            #alignment_length_thres = int(alignment_length_rate_thres * contig_length)
            if evalue <= evalue_thres:
                classified_read_set |= set(read_list)
    if 'root' in classified_read_set:
        classified_read_set.remove('root')
    return classified_read_set
             
# return k longest paths of a graph G including root.
def get_k_longest_paths(G, K):
    assert K != 0
    k_longest_paths = []
    path_nodes = []
    path_nodes = get_path_nodes(G, K)
    end_nodes = get_end_nodes(G)
    for node in end_nodes:
        for weight, path_node in path_nodes[node]:
            add_path(path_nodes, path_node, K, k_longest_paths)
    k_longest_paths.sort()
    return [path for weight, path in k_longest_paths]

def get_end_nodes(G):
    assert G.has_node('root')
    end_nodes = []
    for node in G.nodes_iter():
        if G.out_degree(node) == 0:
            end_nodes.append(node)
    return end_nodes

def test_path_valid(k_longest_paths):
    path_num = len(k_longest_paths)
    for i in range(path_num-1):
        set1 = set(k_longest_paths[i])
        for j in range(i+1, path_num):
            set2 = set(k_longest_paths[j])
            assert not (set1 <= set2 or set1 >= set2)

def add_path_node(path_node, K, path_nodes):
    end = path_node['end']
    path_nodes.setdefault(end, [])
    if len(path_nodes[end]) < K:
        heapq.heappush(path_nodes[end], (path_node['weight'], path_node))
    elif path_node['weight'] > path_nodes[end][0][0]:
        heapq.heapreplace(path_nodes[end], (path_node['weight'], path_node))

def add_path(path_nodes, path_node, K, k_longest_paths):
    if len(k_longest_paths) < K:
        contig = traceback(path_node, path_nodes)
        heapq.heappush(k_longest_paths, (path_node['weight'], contig))
    elif path_node['weight'] > k_longest_paths[0][0]:
        contig = traceback(path_node, path_nodes)
        heapq.heapreplace(k_longest_paths, (path_node['weight'], contig))

def traceback(path_node, path_nodes):
    contig = []
    begin = path_node['begin']
    end = path_node['end']
    kth = path_node['kth']
    while kth != -1:
        contig.append(end)
        current_node = path_nodes[begin][kth][1]
        begin = current_node['begin']
        end = current_node['end']
        kth = current_node['kth']
    contig.append(end)
    contig.reverse()
    return contig

def get_path_nodes(G, K):
    path_nodes = {}
    # add root.
    path_nodes['root'] = []
    path_node = dict(begin='dummy', end='root', weight=0.0, kth=-1)
    path_nodes['root'].append((path_node['weight'], path_node))
    sorted_nodes = nx.topological_sort(G)
    # DP.
    for begin in sorted_nodes:
        for end in G[begin]:
            for i in range(len(path_nodes[begin])):
                path_node = path_nodes[begin][i][1]
                weight = path_node['weight'] + G[begin][end]['weight']
                kth = i
                child_path_node = dict(begin=begin, end=end, weight=weight, kth=kth)
                add_path_node(child_path_node, K, path_nodes)
    return path_nodes

# given mapping and classification result, output sens, fp_rate and ppv.
def output_classification_stat(mapped_read_set, classified_read_set,
                               target_domain, TEST_READ_NUM):
    # TP, FN, FP, TN. 
    counts = get_confusion_mat(mapped_read_set, classified_read_set, TEST_READ_NUM)
    sens = float(counts[0])/(counts[0]+counts[1])  
    fp_rate = float(counts[2])/(counts[2]+counts[3]) 
    if classified_read_set:
        ppv = float(counts[0])/(counts[0]+counts[2])
        out_tuple = (target_domain, counts[0], counts[1], counts[2], counts[3],
                     sens, fp_rate, ppv)
        print '%s %d %d %d %d %.4f %.4e %.4f' % out_tuple
    else:
        ppv = -1
        out_tuple = (target_domain, counts[0], counts[1], counts[2], counts[3],
                     sens, fp_rate, ppv)
        print '%s %d %d %d %d %.4f %.4e %d' % out_tuple


def build_classification_graph(compressed_read_dict, target_domain,
                               fasta_file_name, overlap_thres):
    # target_read_dict is the read dict which will be used to build a graph.
    STATE_THRES = 3
    target_read_num = len(compressed_read_dict)
    G = nx.DiGraph(domain=target_domain)
    target_read_list = compressed_read_dict.values()
    target_read_list.sort(key=lambda read: (read.begin_state, read.end_state))
    for read in target_read_list:
        G.add_node(read.name,
                   begin_state=read.begin_state,
                   end_state=read.end_state,
                   score=read.score,
                   seq = read.seq,
                   members=read.members)
    for i in xrange(target_read_num-1):
        for j in xrange(i+1, target_read_num):
            read1 = target_read_list[i]
            read2 = target_read_list[j]
            read1_length = len(read1.seq)
            read2_length = len(read2.seq)
            pos_overlap_length = \
                get_pos_overlap_length(read1.begin_state, read1.end_state,
                                       read2.begin_state, read2.end_state)
            if pos_overlap_length <= 0:
                continue
            # handle the case that two reads have the same starting state.
            if abs(read1.begin_state - read2.begin_state) <= STATE_THRES:
                seq_overlap_length1 = get_seq_overlap_length(read1.seq, read2.seq)
                seq_overlap_length2 = get_seq_overlap_length(read2.seq, read1.seq)
                if max(seq_overlap_length1, seq_overlap_length2) >= overlap_thres:
                    if seq_overlap_length1 >= seq_overlap_length2:
                        weight = read2.score * \
                                 (read2_length - seq_overlap_length1) / read2_length
                        G.add_edge(read1.name, read2.name,
                                   overlap=seq_overlap_length1, weight=weight)
                    else:
                        weight = read1.score * \
                                 (read1_length - seq_overlap_length1) / read1_length
                        G.add_edge(read2.name, read1.name,
                               overlap=seq_overlap_length2, weight=weight)
                # add an edge between read1 and read2 if they have significant overlap.
                # weight of the edge will be score of read2, which is the child.
            else:
                seq_overlap_length = get_seq_overlap_length(read1.seq, read2.seq)
                if seq_overlap_length >= overlap_thres:
                    weight = read2.score * \
                             (read2_length - seq_overlap_length) / read2_length
                    G.add_edge(read1.name, read2.name,
                               overlap=seq_overlap_length, weight=weight)
    return G

# remove edges that do not help introduce any reads.
def remove_redundant_edges(G):
    for begin_node, end_node in G.edges():
        current_edge_data = G.edge[begin_node][end_node]
        G.remove_edge(begin_node, end_node)
        if not nx.has_path(G, begin_node, end_node):
            G.add_edge(begin_node, end_node, current_edge_data) 

# remove directed cycles in the graph. make it DAG.
def remove_cycles(G):
    while not nx.is_directed_acyclic_graph(G):
        subgraphs = nx.strongly_connected_component_subgraphs(G)
        for subgraph in subgraphs:
            if subgraph.number_of_nodes() > 1:
                edge_index = random.randrange(subgraph.number_of_edges())
                edge = subgraph.edges()[edge_index]
                G.remove_edge(edge[0], edge[1])

# get reads that are not border reads.
def get_target_read_dict(aligned_read_dict, border_read_set):
  target_read_dict = {}
  for read_name in aligned_read_dict:
    if read_name not in border_read_set:
      target_read_dict[read_name] = aligned_read_dict[read_name]
  return target_read_dict

# output classification result.
def output_classified_read_set(classified_read_set, target_domain):
    print '>' + target_domain
    read_set = set()
    for read in classified_read_set:
        read_set.add(read.split('$')[0])	
    for read in sorted(list(read_set)):
        print read

def output_compressed_read_dict(compressed_read_dict):
    for tag in compressed_read_dict:
        print tag, compressed_read_dict[tag].seq

def main(): 
    if len(sys.argv) < 8:
        print >> sys.stderr, 'Usage: <hmmscore file> ' \
                             '<fasta file> <target domain> <overlap thres> ' \
                             '<number of selected paths> <E-value threshold> <working dir>'
        sys.exit(2)
    alignment_file_name = sys.argv[1]
    fasta_file_name = sys.argv[2]
    target_domain = sys.argv[3]
    overlap_thres = int(sys.argv[4])
    N = int(sys.argv[5])
    EVALUE_THRES = float(sys.argv[6])
    working_dir = sys.argv[7]
    aligned_read_dict = get_aligned_read_dict(alignment_file_name, target_domain)
    target_read_dict = aligned_read_dict
    target_read_dict = get_trimmed_aligned_read_dict(aligned_read_dict)
    set_read_seq(fasta_file_name, target_read_dict)
    compressed_read_dict = get_compressed_read_dict(target_read_dict, 'tag')
    if not compressed_read_dict or N == 0:
        print '>' + target_domain
        sys.exit()
    G = build_classification_graph(compressed_read_dict, target_domain, 
                                   fasta_file_name, overlap_thres)
    G.graph['overlap'] = overlap_thres
    remove_redundant_edges(G)
    # remove cycles.
    remove_cycles(G)
    add_root_to_subgraph(G)
    # use number of end nodes if K is invalid.
    if N < 0:
        end_nodes = get_end_nodes(G)
        N = len(end_nodes)
    classifier = classify_top_N_contigs_with_hmmer
    classified_tag_set = classifier(G, N, EVALUE_THRES, working_dir)
    classified_read_set = get_classified_read_set(classified_tag_set, compressed_read_dict)
    output_classified_read_set(classified_read_set, target_domain)

if __name__ == '__main__':
    main()
