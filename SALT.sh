#!/bin/bash
# Copyright (c) 2011 Yuan Zhang, Yanni Sun.
# You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
# pipeline of MetaDomain. 

# input: 
# -m: hmm file;
# -f: fasta file;
# -g: alignment score threshold rate, default: 0.3;
# -k: overlap threshold, default: (average read length) / 2;
# -K: number of contigs generated, default: number of sinks;
# -E: E-value threshold for contigs, default: 1e-6.

# output:
# -o: short read classification in the following output:
# [domain name] 
# [read1]
# [read2]
#.
#.
#.
# [readN1]


usage() {
  echo "./SALT.sh -m <HMMER3 HMM file> -f <fasta file> [options] 
  Options:
    -h:  show this message
    -g:  gamma (alignment score threshold rate, in the range of [0,1], default: 0.3)
    -k:  overlap threshold, default: (average read length) / 2
    -K:  number of contigs generated, default: numbber of sinks
    -E:  E-value threshold for contigs, default: 1e-6
    -o:  output file name, default: stdandard error"
}

hmm=
fasta=
gamma=0.3
k=-1
K=-1
evalue=1e-6
out=
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

while getopts "hm:f:g:k:K:E:o:" OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    m)
      hmm=$OPTARG
      ;;
    f)
      fasta=$OPTARG
      ;;
    g)
      gamma=$OPTARG
      ;;
    k)
      k=$OPTARG
      ;;
    c)
      K=$OPTARG
      ;;
    E)
      evalue=$OPTARG
      ;;
    o)
      out=$OPTARG     
      ;;
    esac
done

if [ "$hmm" == "" ];then
  echo "Please specify the hmm file!"
  usage
  exit
fi

if [ "$fasta" == "" ];then
  echo "Please specify the input fasta file!"
  usage
  exit
fi

if [ `which hmmsearch 2> /dev/null | wc -l` -eq 0 ]; then 
  echo "hmmsearch is not found."; 
  usage
  exit
fi

if [ `$DIR/check_python_packages.py` -eq 1 ];then
  echo "Biopython or NetworkX is not found."
  usage
  exit
fi 

# set k if not given valid value.
if [ $k -le 0 ];then
  k=`awk '/^>/{++num}$1!~/^>/{len+=length}END{print int(len/num/2)}' $fasta`
fi

tmp="salt-tmp"
# generate a temporary folder.
if [ ! -d $tmp ];then
  mkdir $tmp
fi

# generate a list of domains in the input hmm file.
python $DIR/parse_hmm_files.py $hmm $tmp/HMMs
base_fasta=`echo $fasta | awk '{split($1,a,"/"); print a[length(a)]}'`
ls $tmp/HMMs | while read line
do
  hmm_acc=`echo $line | awk '{print substr($1,1,7)}'`
  $DIR/HMMSCORE/hmmscore $tmp/HMMs/$line $fasta $gamma 1 >$tmp/${base_fasta}_${hmm_acc}.hmmscore
  $DIR/analyze_hmmscore_file.py $tmp/${base_fasta}_${hmm_acc}.hmmscore >$tmp/${base_fasta}_${hmm_acc}.K3.hmmscore
  if [ -n "$out" ];then
    $DIR/metadomain.py $tmp/${base_fasta}_${hmm_acc}.K3.hmmscore $fasta ${hmm_acc} \
      $k $K $evalue $tmp >>"$out"
  else 	
    $DIR/metadomain.py $tmp/${base_fasta}_${hmm_acc}.K3.hmmscore $fasta ${hmm_acc} \
      $k $K $evalue $tmp
fi
done
# reset output file if given.
rm -r $tmp
