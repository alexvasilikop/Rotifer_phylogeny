#!/usr/bin/env python3

from sys import argv
from collections import defaultdict

my_fasta     = argv[1]
my_ids	     = argv[2]
my_fasta_out = argv[3]
list_ids = []
ids_seqs=defaultdict(lambda: "")
fasta=defaultdict(lambda: "")

with open (my_ids, "r") as fh_in_ids:

	for line in fh_in_ids:

		line = line.rstrip("\n")

		list_ids.append(">"+line)

id_found=0
id_name=""

with open (my_fasta, "r") as fh_in:

	for line in fh_in:

		line = line.strip()

		if ">" in line:
			id_name = line

		else:
			ids_seqs[id_name]+=line

fh_out=open(my_fasta_out, "w")

for key, seq in ids_seqs.items():

	for ID in list_ids:

		if ID in key:
			fh_out.write(key+"\n"+seq+"\n")

fh_out.close()



