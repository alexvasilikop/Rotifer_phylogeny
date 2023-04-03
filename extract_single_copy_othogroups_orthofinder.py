#!/usr/bin/env python3

import os
import argparse
import re
from collections import defaultdict

# Extract single-copy orthologs from output of orthofinder and print them in new fasta files so that the headers
# have the species name and not the sequence id only (fasta files can then be used directly for alignment and concatenation )
#The script assumes all gene ids are unique (no common ids among proteomes)

def make_out_dir(path):

	'''Check if output dir exists, otherwise make output directory: path_to_out_dir''' 
	if not os.path.isdir(path):
		os.mkdir(path)

def parse_fasta(fileneme, path):

	dict_fasta =defaultdict(lambda: "")
	fh_in_fasta=open(path+fileneme, "r")

	for line in fh_in_fasta:
		line=line.rstrip("\n")

		if line.startswith(">"):
			header=line
		else:
			dict_fasta[header]+=line

	return dict_fasta

def write_fasta(dict_for_fasta, path, OG):

	fh_out=open(path+OG, "w")

	for k,v in dict_for_fasta.items():

		fh_out.write(k+"\n"+v+"\n")

def main():

	parser=argparse.ArgumentParser("Parse output of orthofinder and extract orthogroups with single-copy orthologs into fasta files in which headers have the species name (i.e. for subsequent concatenation).")
	parser.add_argument("orthogroups",                      help="Full path to Orthogroups.tsv")
	parser.add_argument("Single_Copy_Orthologue_Sequences", help="Full path to Single_Copy_Orthologue_Sequences directory")
	parser.add_argument("Output_dir",                       help="Full path to directory to place the output fasta files")

	args=parser.parse_args()

	species_names_ordered=[]
	OG_to_IDs = defaultdict(lambda:[])
	linecounter=0

	#Check if output directory exists and make directory
	make_out_dir(args.Output_dir)

	with open (args.orthogroups, "r") as fh_og:

		for line in fh_og:
			line = line.rstrip("\n")

			if linecounter==0:
				species_names_ordered = line.split("\t")

			else:
				elements = line.split("\t")

				for i in range(1,len(elements)):

					#splitting of appended elements not necessary here as there is one gene copy per OG in single-copy orthologs
					OG_to_IDs[elements[0]].append(elements[i])

			linecounter+=1

	my_single_OG_IDs= os.listdir(args.Single_Copy_Orthologue_Sequences)

	for single_OG_ID in my_single_OG_IDs:
		
		heads_seqs_new = defaultdict(lambda: "")
		heads_seqs=      parse_fasta(single_OG_ID, args.Single_Copy_Orthologue_Sequences)

		og = single_OG_ID.replace(">","").replace(".fa","")
	
		for header in heads_seqs:

			capture_index=1

			for element in OG_to_IDs[og]:

				if header.replace(">","") == element:
					#remove characters fter "." in header including the dot
					heads_seqs_new[">"+re.sub(r"\..+","",species_names_ordered[capture_index])] = heads_seqs[header]
				capture_index+=1

		#writes new fasta in output directory path
		write_fasta(heads_seqs_new, args.Output_dir, single_OG_ID)

##########################################################################
if __name__ == '__main__':
	main()

	



