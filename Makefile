all:
	g++ -O3 -w -o HMMSCORE/hmmscore HMMSCORE/hmmscore.cpp HMMSCORE/viterbi_model.cpp HMMSCORE/DNA_seq.cpp HMMSCORE/hmm_model.cpp
	bash -c "chmod 755 chmod_scripts.sh"
	bash -c "bash chmod_scripts.sh"
