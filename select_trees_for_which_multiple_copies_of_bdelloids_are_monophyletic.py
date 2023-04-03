#!/usr/bin/env python3

import os
from collections import defaultdict
from sys import argv
import re
from ete3 import Tree


# Takes as input trees and corresponding unaligned orthogroup sequences (fasta) for which only one 
#target group of species have multiple copies (e.g. clade Bdelloidea of rotifera); Check that sequences 
#of target clade are monophyletic. Prune trees to keep target sequences and their corresponding fasta files 

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

#######################################################################################
def parse_species_list(file_list_of_species):

	species_clade=defaultdict(lambda: "")

	with open(file_list_of_species, "r") as fh:
		
		for line in fh:
			elements=line.rstrip("\n").split("\t")
			species_clade[elements[0]]=elements[1]

	return species_clade

##################################################################################################
def parse_headers_seqs_target(dict_fasta, ID, clades_dict):

	ids_target=[]

	for i in clades_dict:

		if clades_dict[i] == ID:
			ids_target.append(i)

	head_seqs_target=defaultdict(lambda: "")

	for h,s in dict_fasta.items():

		for i in ids_target:

			if i in h:
				head_seqs_target[(re.sub(r">","", h))]=s

	return head_seqs_target

#############################################################################################################
def main():
	
	file_list_of_species  =argv[1]
	path_to_trees         =argv[2]
	path_to_unaligned_OGs =argv[3]
	group_duplicated      =argv[4]
	outpath               =argv[5]
	treefiles=os.listdir(path_to_trees)
	species_clades=parse_species_list(file_list_of_species)
	
	make_out_dir(outpath)

	for treefile in treefiles:

		print(f"Working on treefile: {path_to_trees+treefile}...")

		t     = Tree(path_to_trees+treefile)
		OG_ID = re.sub(r"_tree.txt", "", treefile)
		heads_seqs_unaligned=parse_fasta(path_to_unaligned_OGs+OG_ID+".fa")
		heads_seqs_target   =parse_headers_seqs_target(heads_seqs_unaligned, group_duplicated, species_clades)

		midpoint_out = t.get_midpoint_outgroup()
		t.set_outgroup(midpoint_out)

		#Check that multiple copies of the target species are monophyletic
		if t.check_monophyly(values=[k for k in heads_seqs_target.keys()], target_attr="name")[1]=='monophyletic':
			
			print("Tree includes monophyletic copies of the target species: "+treefile)

			print("Pruning the tree to remove non-target clade species: "+treefile)

			print("Original tree looks like this:"+"\n")

			print(t)

			ids_to_prune=[x for x in heads_seqs_target.keys()]

			subtree = t.prune(ids_to_prune)

			t.prune(ids_to_prune)
			
			print("Pruned tree"+"\n")
			
			print(t)

			print("Printing subgroup sequences in fasta for tree: "+treefile)

			fh_out_fa = open(outpath+OG_ID+".fa", "w")

			for h,s in heads_seqs_target.items():

				fh_out_fa.write(">"+h+"\n"+s+"\n")

			print("Printing subtree with the sequences of target group for tree: "+treefile)

			t.write(outfile=outpath+OG_ID+".treefile")

			print("Done processing tree: "+treefile)

		else:
			print(f"Seqs of target species in tree {treefile} are not monophyletic")

if __name__ == '__main__':
	main()










