#!/bin/bash
# this file is to install SAT-Assembler in the current folder.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR
chmod 755 hmmer3_pipeline.sh
chmod 755 parse_hmm_files.py 
chmod 755 metadomain.py
chmod 755 check_python_packages.py
chmod 755 analyze_hmmscore_file.py
chmod 755 HMMSCORE/hmmscore 
