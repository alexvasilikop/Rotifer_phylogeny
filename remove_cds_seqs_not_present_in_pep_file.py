#!/usr/bin/env python3

from sys import argv
from collections import defaultdict
import re

#Remove cds (input 2) sequences not present in pep (input 1) file

in_fasta_1=argv[1]
in_fasta_2=argv[2]
out_fasta=argv[3]
heads_seqs_1=defaultdict(lambda:"")
heads_seqs_2=defaultdict(lambda:"")
count_stop_seqs=0

with open(in_fasta_1, "r")as fh_in:

	header=""

	for line in fh_in:

		line=line.rstrip("\n")

		if ">" in line:
			header=line

		else:
			heads_seqs_1[header]+=line

with open(in_fasta_2, "r")as fh_in:

	header=""

	for line in fh_in:

		line=line.rstrip("\n")

		if ">" in line:
			header=line

		else:
			heads_seqs_2[header]+=line


fh_out=open(out_fasta, "w")

for i in heads_seqs_2.keys():

	if i in heads_seqs_1:
		fh_out.write(i+"\n"+heads_seqs_2[i]+"\n")

fh_out.close()






