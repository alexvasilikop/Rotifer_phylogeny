#!/usr/bin/env python3

import os
import argparse
import re
from collections import defaultdict

# Extract single-copy orthologs from output of orthofinder and print them in new fasta files so that the headers correspond to species names
# -> Have the species name and not the sequence id only (fasta files can then be used directly for alignment and concatenation )
# -> Missing copies for some species are allowed (i.e. not all species need to be represented in the orthogroup)
# The script assumes that sequences have unique ids (i.e. no sequence with the same id exists in 2 different proteomes)

def make_out_dir(path):

	'''Check if output dir exists, otherwise make output directory: path_to_out_dir''' 
	if not os.path.isdir(path):
		os.mkdir(path)

def parse_fasta(fileneme, path):

	dict_fasta =defaultdict(lambda: "")
	fh_in_fasta=open(path+fileneme+".fa", "r")

	for line in fh_in_fasta:
		line=line.rstrip("\n")

		if line.startswith(">"):
			header=line
		else:
			dict_fasta[header]+=line

	return dict_fasta

def write_fasta(dict_for_fasta, path, OG):

	fh_out=open(path+OG+".fa", "w")

	for k,v in dict_for_fasta.items():
		fh_out.write(k+"\n"+v+"\n")

def single_copy_check(OG_elements):

	for i in range(1,len(OG_elements)):

		if len(OG_elements[i].split(','))>1:
			return False
	return True


def main():

	parser=argparse.ArgumentParser("Parse output of orthofinder and extract orthogroups with single-copy orthologs into fasta files in which headers have the species name (i.e. for subsequent concatenation).")
	parser.add_argument("Orthogroups",                      help="Full path to Orthogroups.tsv")
	parser.add_argument("Orthogroups_Sequences", 			help="Full path to Orthogroups_Sequences directory")
	parser.add_argument("Output_dir",                       help="Full path to directory to place the output fasta files")
	parser.add_argument("Min_species_number",               help="Min. number of species in Orthogroup (max. 1 gene copy per species)")

	args=parser.parse_args()

	species_names_ordered=[]
	OG_to_IDs = defaultdict(lambda:[])
	linecounter=0

	#Check if output directory exists and make directory
	make_out_dir(args.Output_dir)

	#select only OGs with max; one copy from each species
	with open (args.Orthogroups, "r") as fh_og:

		for line in fh_og:
			line = line.rstrip("\n")

			if linecounter==0:
				species_names_ordered = line.split("\t")

			else:
				elements = line.split("\t")

				#check if OG includes more than one copy in some species
				multiple_copies=single_copy_check(elements)

				if multiple_copies:
					for i in range(1,len(elements)):
						#splitting of appended elements not necessary here as there is one gene copy per OG in single-copy orthologs
						OG_to_IDs[elements[0]].append(elements[i])

			linecounter+=1

	no_single_copy_OGs=0

	for OG_ID in OG_to_IDs:

		heads_seqs_renamed_headers = defaultdict(lambda: "")
		heads_seqs=      parse_fasta(OG_ID, args.Orthogroups_Sequences)

		for header in heads_seqs:
			capture_index=1

			for element in OG_to_IDs[OG_ID]:

				if header.replace(">","") == element:
					#remove characters after "." in header including the dot
					heads_seqs_renamed_headers[">"+species_names_ordered[capture_index]] = heads_seqs[header]
				capture_index+=1


		if len(heads_seqs_renamed_headers.keys())>=int(args.Min_species_number):
			no_single_copy_OGs+=1
			#writes new fasta in output directory path
			write_fasta(heads_seqs_renamed_headers, args.Output_dir, OG_ID)

	print(f"No. of single-copy orthogroups with min. {args.Min_species_number} species: {no_single_copy_OGs}")

##########################################################################
if __name__ == '__main__':
	main()

	



