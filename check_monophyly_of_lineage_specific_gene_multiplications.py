#!/usr/bin/env python3

import os
from collections import defaultdict
from sys import argv
import re
from ete3 import Tree
from statistics import mean


# Check monophyly of lineage-specific gene multiplications (i.e. species-specific multiplications). If monophyletic, 
# Keep only 1 copy according to the following criteria (if first criterion is a tie the next is evaluated)
# 1) Average percent of non-ambiguous coverage for each site (taking into consideration number of residues from other species in the alignment that not ambiguous for this site)
# 2) Longest non-ambiguous sequence
# 3) Select one at random

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
def parse_headers_target(dict_fasta, ID):

	head_seqs_target=defaultdict(lambda: "")

	for h,s in dict_fasta.items():

		if ID in h:
			head_seqs_target[(re.sub(r">","", h))]=s

	return head_seqs_target

#############################################################################################################
def main():
	
	file_list_of_species  =argv[1]
	path_to_trees         =argv[2]
	path_to_alignments    =argv[3]
	path_to_unaligned_OGs =argv[4]
	species_duplicated    =argv[5]
	outpath               =argv[6]
	treefiles=os.listdir(path_to_trees)
	species_clades=parse_species_list(file_list_of_species)
	
	make_out_dir(outpath)

	for treefile in treefiles:

		print(f"Working on treefile: {path_to_trees+treefile}...")

		t     = Tree(path_to_trees+treefile)
		OG_ID = re.sub(r"_tree.txt", "", treefile)
		heads_seqs_aligned  =parse_fasta(path_to_alignments+OG_ID+".fa")
		heads_seqs_unaligned=parse_fasta(path_to_unaligned_OGs+OG_ID+".fa")
		heads_seqs_target   =parse_headers_target(heads_seqs_aligned, species_duplicated)

		midpoint_out = t.get_midpoint_outgroup()
		t.set_outgroup(midpoint_out)

		#Check that multiple copies of the target species are monophyletic
		if t.check_monophyly(values=[k for k in heads_seqs_target.keys()], target_attr="name")[1]=='monophyletic':
			print(heads_seqs_target.keys())
			print("Tree includes monophyletic copies of the target species: "+treefile)

			#Select copy with the highest average non-ambiguous coverage of residues in the alignment
			IDs_scores=defaultdict(lambda: 0)
			counter=0

			for h in heads_seqs_target.keys():
				print(f"Calculating average non-ambiguous coverage of copy: "+h)

				index=0
				percent_covered_all_sites=[]

				for r in heads_seqs_target[h]:

					total_covered=0

					if r!="X" and r!="-":

						for h2 in heads_seqs_aligned.keys():

							if (h2 != h) and (species_duplicated not in h2):

								if heads_seqs_aligned[h2][index]!="X" and heads_seqs_aligned[h2][index]!="-":
									total_covered+=1

					percent_covered_all_sites.append(total_covered/(len(heads_seqs_aligned)-len(heads_seqs_target)))
					index+=1

				IDs_scores[h]=mean(percent_covered_all_sites)
				print("Number of ambiguous characters for copy: {ambiguous}".format(ambiguous=len(re.findall(r"-|X", heads_seqs_aligned[h]))))
			print(IDs_scores)

			score_max=0
			ID_max=""

			if len(set(IDs_scores.values())) < len(IDs_scores.keys()):
				print("Found equal coverage scores -> Selecting best copy based on non-ambiguous length")

				#return keys with equal average coverage
				max_list_keys=[k for k, v in IDs_scores.items() if v == max(IDs_scores.values())]
				print(max_list_keys)

				if len(max_list_keys)==1:
					score_max=max(IDs_scores.values())
					ID_max=max_list_keys[0]

				else:
					max_non_ambiguous_len=0

					#if both are of equal non-ambiguous length the last one is selected at random
					for ID in max_list_keys:
						if ((len(heads_seqs_aligned[ID])-len(re.findall(r"-|X", heads_seqs_aligned[ID])))) > max_non_ambiguous_len:
							score_max=(len(heads_seqs_aligned[ID])-len(re.findall(r"-|X", heads_seqs_aligned[ID])))
							ID_max=ID

			else:
				for ID in IDs_scores.keys():

					if IDs_scores[ID]>score_max:
						score_max=IDs_scores[ID]
						ID_max=ID

			fh_out=open(outpath+OG_ID+".fa", "w")

			for header,sequence in heads_seqs_unaligned.items():

				if species_duplicated in header:

					if ID_max in header:
						fh_out.write(">"+header+"\n"+sequence+"\n")

				else:
					fh_out.write(">"+header+"\n"+sequence+"\n")

			fh_out.close()

		else:
			print(f"Seqs of target species in tree {treefile} are not monophyletic")

if __name__ == '__main__':
	main()










