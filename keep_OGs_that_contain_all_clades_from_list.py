#!/usr/bin/env python3

import os
from collections import defaultdict
import shutil
import glob
from sys import argv
import re

def make_out_dir(path):

	'''Check if output dir exists, otherwise make output directory: path_to_out_dir''' 
	if not os.path.isdir(path):
		os.mkdir(path)

def parse_fasta(fileneme):

	dict_fasta =defaultdict(lambda: "")
	fh_in_fasta=open(fileneme, "r")

	for line in fh_in_fasta:
		line=line.rstrip("\n")

		if line.startswith(">"):
			header=line
		else:
			dict_fasta[header]+=line

	return dict_fasta

def parse_species_list(file_list_of_species):

	species_clade=defaultdict(lambda: "")

	with open(file_list_of_species, "r") as fh:
		
		for line in fh:
			elements=line.rstrip("\n").split("\t")
			print(elements)
			species_clade[elements[0]]=elements[1]

	return species_clade

def main():
	
	file_list_of_species=argv[1]
	outpath=argv[2]
	files=glob.glob("*.fa")
	species_clades=parse_species_list(file_list_of_species)
	clades=set(species_clades.values())
	print(f"Clades: {clades}")
	print(f"No. Clades: {len(clades)}")
	
	make_out_dir(outpath)

	for f in files:
		print(f)
		no_clades_in_file=defaultdict(lambda: 0)

		heads_seqs=parse_fasta(f)

		for h in heads_seqs.keys():
			if re.sub(">","",h) in species_clades.keys():
				no_clades_in_file[species_clades[re.sub(">","",h)]]+=1

		print(f"File: {f}")
		print(f"Clades: {no_clades_in_file.keys()}")
		print(f"No. of clades: {len(no_clades_in_file)}")

		if len(no_clades_in_file)==len(clades):

			shutil.copy(f,outpath)

if __name__ == '__main__':
	main()










