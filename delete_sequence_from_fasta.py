#!/usr/bin/env python3

#Delete a specific sequence from fasta file
#Help message: python3 delete_sequence_from_fasta -h

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
def filter_fasta(dict_fasta,species_name):
	
	try:
		del dict_fasta[species_name]
	except KeyError:
		print(species_name+" not found in fasta!")
	finally:
		return dict_fasta

##################################################
def main():

	parser=argparse.ArgumentParser(description="Remove specific sequence from fasta. If species not in fasta the output is the same unfiltered fasta")
	parser.add_argument("in_fasta", help="Input fasta file")
	parser.add_argument("out_fasta", help="Output fasta file (filtered)")
	parser.add_argument("species", help="Species to remove from fasta")

	args=parser.parse_args()
	seqs_headers=parse_fasta(args.in_fasta)
	seqs_headers_new=filter_fasta(seqs_headers,args.species)
	write_fasta(seqs_headers_new, args.out_fasta)

###################################################################
if __name__=='__main__':
	main()