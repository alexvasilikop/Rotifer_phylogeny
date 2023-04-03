#!/usr/bin/env python3

#Delete all sequences matching a specified string from fasta file
#Help message: python3 delete_sequence_with_pattern_from_fasta -h

import argparse
from collections import defaultdict
import re

#################################################
def parse_fasta(filename):

	seqs_heads=defaultdict(lambda:"")
	header=""
	with open (filename, "r") as fh:

		for line in fh:
			line=line.rstrip("\n")

			if ">" in line:
				header=line
			else:
				seqs_heads[re.sub(">","",header)]+=line

	return seqs_heads

#############################################################
def write_fasta(dict_fasta, filename):

	fh_out=open(filename, "w")

	for i,j in dict_fasta.items():

		fh_out.write(">"+i+"\n"+j+"\n")

###########################################################
def filter_fasta(dict_fasta,pattern):
	
	filtered_dict=defaultdict(lambda:"")

	for key,seq in dict_fasta.items():
		if pattern not in key:
			filtered_dict[key]=seq

	return filtered_dict

##################################################
def main():

	parser=argparse.ArgumentParser(description="Remove all sequences matching a specified string from fasta")
	parser.add_argument("in_fasta", help="Input fasta file")
	parser.add_argument("out_fasta", help="Output fasta file (filtered)")
	parser.add_argument("string", help="String pattern required for a sequence to be deleted from the fasta")

	args=parser.parse_args()
	seqs_headers=parse_fasta(args.in_fasta)
	seqs_headers_new=filter_fasta(seqs_headers,args.string)
	write_fasta(seqs_headers_new, args.out_fasta)

###################################################################
if __name__=='__main__':
	main()