#!/usr/bin/env python3

#Filter trinity output to keep the longest transcript
#Usage: python3 keep_longest_trinity_transcript.py in_assembly out_assembly

from sys import argv
from collections import defaultdict
import re

in_assembly=argv[1]
out_assembly=argv[2]

def empty_dict():
	return ""

header=""
head_seqs=defaultdict(empty_dict)
head_seqs_new=defaultdict(empty_dict)

with open (in_assembly, "r") as fh_in:

	for line in fh_in:

		line=line.rstrip("\n")

		if line.startswith(">"):

			elements = line.split(" ")
			gene_is_info=elements[0].split("_")
			header = gene_is_info[0]+"_"+gene_is_info[1]+"_"+gene_is_info[2]+"_"+gene_is_info[3]+" "+gene_is_info[4]+" "+elements[1]
			print(header)

		else:
			head_seqs[header]+=line


for h,s in head_seqs.items():

	elements=h.split(" ")

	if elements[0] in head_seqs_new.keys():

		if len(s) > len(head_seqs_new[elements[0]]):

			head_seqs_new[elements[0]]=s

	else:
		head_seqs_new[elements[0]]=s

fh_out=open(out_assembly, "w")
for filtered_h, filtered_s in head_seqs_new.items():

	fh_out.write(filtered_h+"\n"+filtered_s+"\n")




	







