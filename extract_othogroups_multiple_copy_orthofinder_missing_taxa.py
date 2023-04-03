#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict
import glob
import os

# Extract orthogroups from output of orthofinder (maximum 2 copies per species) and print them in new fasta files so that the headers correspond to species names
# -> Have also the species name and not the sequence id only (fasta files can then be used directly for alignment and concatenation )
# -> Missing copies for some species are allowed (i.e. not all species need to be represented in the orthogroup)
# The script assumes that sequences have unique ids (i.e. no sequence with the same id exists in 2 different proteomes)

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

	fh_out=open(path+OG+".fa", "w")

	for k,v in dict_for_fasta.items():
		fh_out.write(k+"\n"+v+"\n")

def check_max_no_copies_per_species(dict_seqs, max_copies):

	gene_copy_counts=defaultdict(lambda: 0)

	for h in dict_seqs.keys():
		elements_h=h.split(" ")
		gene_copy_counts[elements_h[0]]+=1

	counts = gene_copy_counts.values()
	print(counts)

	for count in counts:
		if count > int(max_copies):
			return False
	return True

def main():

	parser=argparse.ArgumentParser("Parse output of orthofinder and extract orthogroups into fasta files in which headers also have the species name (missing taxa in some orthogroups are allowed).")
	parser.add_argument("Orthogroups",                      help="Full path to Orthogroups.tsv")
	parser.add_argument("Orthogroups_Sequences", 			help="Full path to Orthogroups_Sequences directory")
	parser.add_argument("Output_dir",                       help="Full path to directory to place the output fasta files")
	parser.add_argument("Max_gene_copies",   	            help="Max no. of gene copies per species in orthogroup")

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

				for i in range(1,len(elements)):
					#splitting of appended elements not necessary here as there is one gene copy per OG in single-copy orthologs
					OG_to_IDs[elements[0]].append(elements[i])

			linecounter+=1

	my_input_OG = glob.glob("*.fa")

	for OG_ID in OG_to_IDs:

		heads_seqs_renamed_headers = defaultdict(lambda: "")
		heads_seqs=      parse_fasta(OG_ID+".fa", args.Orthogroups_Sequences)

		for header in heads_seqs:

			capture_index=1

			for element in OG_to_IDs[re.sub(".fa","",OG_ID)]:

				if header.replace(">","") in element:
					#remove characters after "." in header including the dot
					heads_seqs_renamed_headers[">"+species_names_ordered[capture_index]+" "+header.replace(">","")] = heads_seqs[header]
				capture_index+=1

		if check_max_no_copies_per_species(heads_seqs_renamed_headers, args.Max_gene_copies):

			print(OG_ID)
			if len(heads_seqs_renamed_headers)>=2:

				#writes new fasta in output directory path
				write_fasta(heads_seqs_renamed_headers, args.Output_dir, OG_ID)


##########################################################################
if __name__ == '__main__':
	main()

	



