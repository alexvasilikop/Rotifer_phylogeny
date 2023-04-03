#!/usr/bin/env python3

import os
from collections import defaultdict
import shutil
import glob
from sys import argv
import re

#Keep OGs that have exactly one copy in all species outside the target clade but have multiple copies in 
#the species of the target clade (at least 2 species from the target group with at least 2 copies each must be present)

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
	files=glob.glob("*.fa")
	species_clades=parse_species_list(file_list_of_species)
	clades=set(species_clades.values())
	print(f"Clades: {clades}")
	print(f"No. Clades: {len(clades)}")
	print(f"Finding orthogroups for which only {ID_multi_copy} group is multi-copy and the rest of groups single-copy")
	
	make_out_dir(outpath)

	for f in files:
		no_clades_in_file=defaultdict(lambda: 0)

		heads_seqs=parse_fasta(f)
		no_species_target_group           =defaultdict(lambda: 0)
		no_copies_per_species_target_group=defaultdict(lambda: 0)
		no_species_other_groups           =defaultdict(lambda: 0)
		no_copies_per_species_other_groups=defaultdict(lambda: 0)

		for h in heads_seqs.keys():

			ID=re.sub(r">","",h)
			ID2=re.sub(r" .*","",ID)
			if ID2 in species_clades.keys():

				if species_clades[ID2] == ID_multi_copy:
					no_species_target_group[ID2]+=1
					no_copies_per_species_target_group[ID2]+=1

				else:
					no_copies_per_species_other_groups[ID2]+=1


		#Check that target group is represented by more than 1 species
		if len(no_species_target_group) > 1:

			no_multi=0
			for i in no_copies_per_species_target_group.values():
				if i>1:
					no_multi+=1

			#Check that at least two species in the target group have at least 2 copies
			if no_multi >1:

				no_multi_non_target=0
				for i in no_copies_per_species_other_groups.values():
					if i>1:
						no_multi_non_target+=1

				#Check that species not belonging to target group do not have more than 1 copy
				if no_multi_non_target==0:

					print(f)
					print(f"Number of target group species: {len(no_species_target_group)}")
					print(f"Number of copies in target group species: {no_copies_per_species_target_group.values()}")
					print(f"Number of copies in non-target group species: {no_copies_per_species_other_groups.values()}")
						
					shutil.copy(f,outpath)


if __name__ == '__main__':
	main()










