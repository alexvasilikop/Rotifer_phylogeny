#!/usr/bin/env python3

import os
from collections import defaultdict
from sys import argv
import re
from ete3 import Tree

#Select the best subgroup (single-copy orthogroup) for these OGs that the target clade (e.g. Bdelloidea) has more than one subgroup and
#print combined fasta to a new directory

def make_out_dir(path):

	'''Check if output dir exists, otherwise make output directory: path_to_out_dir''' 
	if not os.path.isdir(path):
		os.mkdir(path)

###########################################################################################""
def parse_fasta(filename):

	dict_fasta =defaultdict(lambda: "")
	fh_in_fasta=open(filename, "r")

	for line in fh_in_fasta:
		line=line.rstrip("\n")

		if line.startswith(">"):
			header=re.sub(r">","",line)
		else:
			dict_fasta[header]+=line

	return dict_fasta

##################################################################################################################
def parse_species_list(file_list_of_species):

	species_clade=defaultdict(lambda: "")

	with open(file_list_of_species, "r") as fh:
		
		for line in fh:
			elements=line.rstrip("\n").split("\t")
			species_clade[elements[0]]=elements[1]

	return species_clade

#############################################################################################################
def main():
	
	species_list          =argv[1]
	path_to_trees         =argv[2]
	path_to_unaligned_OGs_bdelloid_subgroups=argv[3]
	path_to_unaligned_OGs =argv[4]
	outpath               =argv[5]
	target_clade          =argv[6]
	treefiles=os.listdir(path_to_trees)
	original_fasta=os.listdir(path_to_unaligned_OGs)
	best_seqs_per_OG=defaultdict(lambda: {})
	heads_seqs_all_OG=defaultdict(lambda: {})
	subtrees_lengths=defaultdict(lambda: 0)
	
	make_out_dir(outpath)

	#parse species-clades list
	species_clades=parse_species_list(species_list)

	#store original fasta sequences for each ID
	for i in original_fasta:

		ID=re.sub(r".fa", "", i)

		heads_seqs_all_OG[ID]=parse_fasta(path_to_unaligned_OGs+i)

	#Select one subgroup of clade for each OG
	for treefile in treefiles:

		#Get OG ID
		OG_ID=re.sub(r".fa.orthosnap...tre", "", treefile)
		orthosnap_ID1=re.sub(rf"{OG_ID}.", "", treefile)
		orthosnap_ID=re.sub(r".tre", "", orthosnap_ID1)
		print(OG_ID)
		print(orthosnap_ID)

		print(f"Working on treefile: {path_to_trees+treefile}...")

		t = Tree(path_to_trees+treefile, format=5)

		#Read fasta subgroup
		try:
			heads_seqs_unaligned_bdelloids=parse_fasta(path_to_unaligned_OGs_bdelloid_subgroups+OG_ID+"."+orthosnap_ID+".fa")

		except:
			print("Error in opening file: "+path_to_unaligned_OGs_bdelloid_subgroups+OG_ID+"."+orthosnap_ID+".fa")
			print("Double check paths and files!")
		
		no_species_for_OG=0
		current_total_tree_length=0

		#count number of species in tree
		for n in t.traverse():

			if not n.is_root():
				current_total_tree_length+=n.dist

			if n.is_leaf():
				no_species_for_OG+=1

		print(f"Treefile has the following no. of species: {no_species_for_OG}")
		print(f"Treefile has total tree length: {current_total_tree_length}")

		if OG_ID in best_seqs_per_OG.keys():
		
			#Update if new subgroup has more species
			if len(best_seqs_per_OG[OG_ID]) < no_species_for_OG:

				best_seqs_per_OG[OG_ID]=heads_seqs_unaligned_bdelloids

			#Update if new subgroup has shorter total tree length
			elif subtrees_lengths[OG_ID] > current_total_tree_length:

				subtrees_lengths[OG_ID] = current_total_tree_length
				best_seqs_per_OG[OG_ID]=heads_seqs_unaligned_bdelloids

		else:
			#store bdelloid sequences in subgroup to global best dictionary
			best_seqs_per_OG[OG_ID]=heads_seqs_unaligned_bdelloids


	for ID in heads_seqs_all_OG:

		print("Printing final sequences in fasta for orthogroup: "+ID)

		fh_out_fa=open(outpath+ID+".fa", "w")

		for h in heads_seqs_all_OG[ID]:

			species=re.sub(r"\|.*","", h)
			clade=species_clades[species]

			if clade != target_clade:
				fh_out_fa.write(">"+h+"\n"+heads_seqs_all_OG[ID][h]+"\n")

			else:
				if h in best_seqs_per_OG[ID]:
					fh_out_fa.write(">"+h+"\n"+heads_seqs_all_OG[ID][h]+"\n")

		fh_out_fa.close()
		print("Done processing orthogroup: "+ID)

if __name__ == '__main__':
	main()










