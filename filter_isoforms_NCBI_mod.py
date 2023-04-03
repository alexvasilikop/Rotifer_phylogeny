#!/usr/bin/env python3 

from sys import argv
from collections import defaultdict
import re

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
				header = line

			else:
				head_seqs_cds[header]+=line

with open(pep_unfiltered, "r") as fh_pep:

	header = ""
	for line in fh_pep:

			line = line.rstrip("\n")

			if line.startswith(">"):
				elements=line.split(" ")
				header = ">"+re.sub(r"gene=", "", elements[2])
			else:
				head_seqs_pep[header]+=line

head_seqs_cds_filtered=defaultdict(lambda: "")
head_seqs_pep_filtered=defaultdict(lambda: "")

#Filter CDS and peptides by length and assign locus ID as header
for h,s in head_seqs_cds.items():
	
		elements=h.split(" ")
		l_id=">"+re.sub(r"gene=", "", elements[1])

		if l_id in head_seqs_cds_filtered.keys():

			if len(s)>=len(head_seqs_cds_filtered[l_id]):
				head_seqs_cds_filtered[l_id]=s
				head_seqs_pep_filtered[l_id]=head_seqs_pep[l_id]
		else:
			head_seqs_cds_filtered[l_id]=s
			head_seqs_pep_filtered[l_id]=head_seqs_pep[l_id]

#Print filtered files
fh_cds_out=open(cds_filtered, "w")
for i,j in head_seqs_cds_filtered.items():
	fh_cds_out.write(i+"\n"+j+"\n")
fh_cds_out.close()

fh_pep_out=open(pep_filtered, "w")
for i,j in head_seqs_pep_filtered.items():
	fh_pep_out.write(i+"\n"+j+"\n")
fh_pep_out.close()

