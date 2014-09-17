#!/usr/bin/python
import sys
from pprint import pprint

def get_read_domain_dict(in_file_name):
    read_domain_dict = {}
    with open(in_file_name, 'Ur') as f:
        for line in f:
            if line.strip():
                items = line.rstrip().split()
                read_name = items[0][1:]
                domain_name = items[1].split('=')[1][0:7]
                score = float(items[2].split('=')[1])
                alignment_aa_length = int(items[6].split('=')[1])
                read_domain_dict.setdefault(read_name, [])
                read_domain_dict[read_name].append((score, domain_name))
    return read_domain_dict

def get_average_read_domain_num(read_domain_dict):
    read_num = len(read_domain_dict)
    read_domain_num = 0
    for read_name in read_domain_dict:
        read_domain_num += len(read_domain_dict[read_name])
    print read_num, read_domain_num, float(read_domain_num) / read_num

def output_read_domain_num(read_domain_dict):
    for read_name in read_domain_dict:
        print read_name, len(read_domain_dict[read_name])
   
def get_top_domains_for_reads(read_domain_dict, K):
    read_top_domain_dict = {}
    for read_name in read_domain_dict:
        read_domain_dict[read_name].sort()
        read_top_domain_dict[read_name] = \
            set([domain for score, domain in read_domain_dict[read_name][-K:]])
    return read_top_domain_dict

def output_new_hmmscore_file(in_file_name, read_top_domain_dict):
    with open(in_file_name, 'Ur') as f:
        for line in f:
            if line.strip():
                items = line.rstrip().split()
                read_name = items[0][1:]
                domain_name = items[1].split('=')[1][0:7]
                if domain_name in read_top_domain_dict[read_name]:
                   print line.rstrip() 

def main():
    if len(sys.argv) < 2:
        print >> sys.stderr, 'Usage: <hmmscore file>'
        sys.exit(2)
    read_domain_dict = get_read_domain_dict(sys.argv[1]) 
    read_top_domain_dict = get_top_domains_for_reads(read_domain_dict, 3);
    output_new_hmmscore_file(sys.argv[1], read_top_domain_dict)

if __name__ == '__main__':
    main()
    

