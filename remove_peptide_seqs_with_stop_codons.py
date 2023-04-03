#!/usr/bin/env python3

from sys import argv
from collections import defaultdict
import re

#Remove peptide sequences with stop codons ("stop is identified as . character")

in_fasta=argv[1]
out_fasta=argv[2]
heads_seqs=defaultdict(lambda:"")
count_stop_seqs=0

with open(in_fasta, "r")as fh_in:

	header=""

	for line in fh_in:

		line=line.rstrip("\n")

		if ">" in line:
			header=line

		else:
			heads_seqs[header]+=line


fh_out=open(out_fasta, "w")

for i,j in heads_seqs.items():

	if "." in j:
		count_stop_seqs+=1

	else:
		fh_out.write(i+"\n"+j+"\n")

fh_out.close()






