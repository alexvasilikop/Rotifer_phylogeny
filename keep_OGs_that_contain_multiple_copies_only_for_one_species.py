#!/usr/bin/env python3

import os
from collections import defaultdict
import shutil
import glob
from sys import argv
import re

#Keep OGs that have multiple copies in only one species but are single-copy in all other species

def make_out_dir(path):

	'''Check if output dir exists, otherwise make output directory: path_to_out_dir''' 
	if not os.path.isdir(path):
		os.mkdir(path)

def parse_fasta(filename):

	dict_fasta =defaultdict(lambda: "")
	fh_in_fasta=open(filename, "r")

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
			species_clade[elements[0]]=elements[1]

	return species_clade

def main():
	
	file_list_of_species=argv[1]
	outpath=argv[2]
	ID_multi_copy=argv[3]
	max_no_copies=argv[4]
	files=glob.glob("*.fa")
	species_clades=parse_species_list(file_list_of_species)
	clades=set(species_clades.values())
	print(f"Clades: {clades}")
	print(f"No. Clades: {len(clades)}")
	print(f"Finding orthogroups for which only species {ID_multi_copy} group is multi-copy and the rest of species single-copy")
	
	make_out_dir(outpath)

	for f in files:
		no_clades_in_file=defaultdict(lambda: 0)

		heads_seqs=parse_fasta(f)
		no_copies_species_target=0
		no_copies_per_species=defaultdict(lambda: 0)

		for h in heads_seqs.keys():

			ID=re.sub(r">","",h)
			ID2=re.sub(r" .*","",ID)
			if ID2 in species_clades.keys():

				if ID2 == ID_multi_copy:
					no_copies_species_target+=1

				else:
					no_copies_per_species[ID2]+=1

		#Check that target species is represented by more than 1 copy (up to a maximum)
		if no_copies_species_target > 1 and no_copies_species_target <= int(max_no_copies):

				no_multi_non_target=0
				for i in no_copies_per_species.values():
					if i>1:
						no_multi_non_target+=1

				#Check that other species do not have more than 1 copy
				if no_multi_non_target==0:

					print(f)
					print(f"Target species copy number: {no_copies_species_target}")
					print(f"Number of copies in non-target species: {no_copies_per_species.values()}")
						
					shutil.copy(f,outpath)
#####################################################################

if __name__ == '__main__':
	main()










