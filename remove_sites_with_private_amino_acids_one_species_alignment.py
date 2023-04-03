#!/usr/bin/env python3

import argparse
from collections import defaultdict
import re

#################################################
def read_fasta(filename):

	#parse fasta in dictionary
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

	#write output fasta file
	fh_out=open(filename, "w")

	for i,j in dict_fasta.items():

		fh_out.write(">"+i+"\n"+j+"\n")

############################################################################
def main():

	parser=argparse.ArgumentParser(description="Remove sites that have unique amino acids for a specified species (i.e. unique/private amino-acid or nuc residue for this species) from a multiple sequence alignment.")
	parser.add_argument("Alignment_fasta", help="Alignment in fasta format")
	parser.add_argument("Species_name", help="Species name in the fasta file")
	parser.add_argument("Output_fasta_name", help="Name of filtered alignment file in fasta format")
	parser.add_argument("Data_type", help="nuc (for nucleotides) or aa (for amino acids)")
	parser.add_argument("Other_singletons", help="If there is a unique amino acid in a column for the specified species, specify whether to keep the column or not depending on if there are singletons for other species in the same column (yes/no)")
	args=parser.parse_args()

	sequences_and_headers=read_fasta(args.Alignment_fasta)
	sequences_and_headers_new=defaultdict(lambda:"")
	no_singleton_sites=0
	total_sites=0
	ambiguous=""

	if args.Data_type == "aa":
		ambiguous="X?-"
	elif args.Data_type == "nuc":
		ambiguous="N?-"
	else:
		print("Unknown data type please specify \"aa\" or \"nuc\"")
		exit(1)

	#check input species name in fasta
	if args.Species_name not in sequences_and_headers:
		print("Species not present in Fasta file!")
		exit(1)

	for index in range(0,len(sequences_and_headers[args.Species_name])):

		residue=sequences_and_headers[args.Species_name][index]

		if residue in ambiguous:
			for species in sequences_and_headers.keys():
					sequences_and_headers_new[species]+=sequences_and_headers[species][index]
			total_sites+=1

		else:
			shared=0
			singleton_others=0

			for species in sequences_and_headers.keys():
				if species != args.Species_name:
					if sequences_and_headers[species][index] == residue:
						shared+=1

			if shared>0:
				for species in sequences_and_headers.keys():
					sequences_and_headers_new[species]+=sequences_and_headers[species][index]
				total_sites+=1

			else:
				#count number of other residues on the same column
				no_residues_others=defaultdict(lambda: 0)
				for species in sequences_and_headers.keys():
					if species != args.Species_name:
						if sequences_and_headers[species][index] not in ambiguous:
							no_residues_others[sequences_and_headers[species][index]]+=1

				#add columns to output if singletons are allowed for other species in the same column otherwise discard columns
				if args.Other_singletons == "no":
						no_singleton_sites+=1

				elif args.Other_singletons == "yes":
					#remove sites for which only Seison is singleton and keep others
					if 1 not in no_residues_others.values():
						no_singleton_sites+=1

					else:
						for species in sequences_and_headers.keys():
							sequences_and_headers_new[species]+=sequences_and_headers[species][index]
						total_sites+=1
						
				else:
					print("Unspecified option for unique residues in other species of the same column: please specify: yes or no")
					exit(1)

	write_fasta(sequences_and_headers_new, args.Output_fasta_name)
	print("No. removed sites: "+str(no_singleton_sites))
	print("No. remaining sites: "+str(total_sites))

##################################################
if __name__ == '__main__':

	main()
	


