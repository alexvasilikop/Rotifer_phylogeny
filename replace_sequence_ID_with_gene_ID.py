#!/usr/bin/env python3

from sys import argv
from collections import defaultdict
import re

in_fasta=argv[1]
in_translation_table=argv[2]
out_fasta=argv[3]

mRNA_to_gene=defaultdict(lambda: "")

with open(in_translation_table, "r") as fh_table:

	for line in fh_table:

		line = line.rstrip("\n")
		elements = line.split(" ")
		mRNA_to_gene[elements[1]]=elements[0]


seqs_headers=defaultdict(lambda: "")
header=""

with open(in_fasta, "r") as fh_fasta_in:

	for line in fh_fasta_in:

		line = line.rstrip("\n")

		if line.startswith(">"):

			header=">"+mRNA_to_gene[re.sub(">","",line)]

		else:
			seqs_headers[header]+=line


fh_out=open(out_fasta, "w")

for i,j in seqs_headers.items():

	fh_out.write(i+"\n"+j+"\n")






