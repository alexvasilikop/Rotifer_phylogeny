#!/usr/bin/env python3 

from sys import argv
from collections import defaultdict
import re

#Filter CDS and PROTEIN files downloaded from NCBI based on locus_tag -> keep one sequence per locus tag
#If multiple cds and corresponding proteins for a locus then the longest cds (and corresponding peptide) are selected
#Also filter out cds and protein predictions corresponding to pseusogenes
#Usage python3 filter_isoforms_NCBI.py in_cds.fa in_pep.fa out_cds.fa out_pep.fa

cds_unfiltered = argv[1]
pep_unfiltered = argv[2]
cds_filtered = argv[3]
pep_filtered = argv[4]
head_seqs_cds = defaultdict(lambda: "")
head_seqs_pep = defaultdict(lambda: "")

with open(cds_unfiltered, "r") as fh_cds:

	header =""
	pseudo=0
	for line in fh_cds:
			line = line.rstrip("\n")

			if line.startswith(">"):

				if "pseudo=true" in line:
					pseudo=1
					next
				else:
					pseudo=0
				header = re.sub(r"[\[\]]", "", line)
				header = re.sub(r"\|", "_", header)
				elements=header.split(" ")

				for i in elements:
					if "protein_id" in i:
						id_found=i
						header= " ".join(elements[:2])
						header+=" "+id_found
			else:
				if pseudo==0:
					head_seqs_cds[header]+=line

with open(pep_unfiltered, "r") as fh_pep:

	header = ""
	for line in fh_pep:

			line = line.rstrip("\n")

			if line.startswith(">"):
				elements=line.split(" ")
				header = elements[0]
			else:
				head_seqs_pep[header]+=line

head_seqs_cds_filtered=defaultdict(lambda: "")
head_seqs_pep_filtered=defaultdict(lambda: "")

#Filter CDS and peptides by length and assign locus ID as header
for h,s in head_seqs_cds.items():
	
	if "protein_id=" in h:
		elements=h.split(" ")
		l_id=">"+re.sub(r"locus_tag=", "", elements[1])
		p_id=re.sub(r"protein_id=", "", elements[2])
		header_p=">"+p_id

		if l_id in head_seqs_cds_filtered.keys():

			if len(s)>=len(head_seqs_cds_filtered[l_id]):
				head_seqs_cds_filtered[l_id]=s
				head_seqs_pep_filtered[l_id]=head_seqs_pep[header_p]
		else:
			head_seqs_cds_filtered[l_id]=s
			head_seqs_pep_filtered[l_id]=head_seqs_pep[header_p]

#Print filtered files
fh_cds_out=open(cds_filtered, "w")
for i,j in head_seqs_cds_filtered.items():
	fh_cds_out.write(i+"\n"+j+"\n")
fh_cds_out.close()

fh_pep_out=open(pep_filtered, "w")
for i,j in head_seqs_pep_filtered.items():
	fh_pep_out.write(i+"\n"+j+"\n")
fh_pep_out.close()

