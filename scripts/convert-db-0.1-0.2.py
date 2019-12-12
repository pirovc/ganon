#!/usr/bin/env python3

import os, sys, pickle, gzip, shutil

if len(sys.argv)==1:
	print("Conversion script for ganon database files from version 0.1.X to 0.2.X")
	print("Usage: ./convert-db-0.1-0.2.py db-prefix-in [db-prefix-out]")
	print("If no db-prefix-out is provided, the input database will be overwritten")
	sys.exit()

db_prefix_in = sys.argv[1]
overwrite=False
if len(sys.argv)==3:
	db_prefix_out = sys.argv[2]
else:
	db_prefix_out = db_prefix_in
	overwrite=True

file_nodes_in = db_prefix_in+".nodes"
file_bins_in = db_prefix_in+".bins"
file_map_in = db_prefix_in+".map"
file_filter_in = db_prefix_in+".filter"

file_tax_out = db_prefix_out+".tax"
file_gnn_out = db_prefix_out+".gnn"
file_map_out = db_prefix_out+".map"
file_ibf_out = db_prefix_out+".ibf"


info = pickle.load(open(file_nodes_in, 'rb'))

# write .tax
tax = open(file_tax_out, "w")
print('\n'.join(['\t'.join([node,tax[0],tax[2],tax[1]]) for node,tax in info['filtered_nodes'].items()]), file=tax)
tax.close()

# generate .gnn
bins = []
with open(file_bins_in, 'r') as file:
    for line in file:
        bins.append(line.strip())

gnn = {'kmer_size': info['kmer_size'], 
'hash_functions': info['hash_functions'], 
'number_of_bins': info['number_of_bins'], 
'rank': info['rank'],
'bin_length': info['bin_length'],
'fragment_length': info['fragment_length'],
'overlap_length': info['overlap_length'],
'bins': bins }
with gzip.open(file_gnn_out, 'wb') as f: pickle.dump(gnn, f)

# make .map uniq
ganon_map = {}
with open(file_map_in, 'r') as file:
    for line in file:
       target, binid = line.rstrip().split("\t")
       ganon_map[binid] = target

with open(file_map_out, 'w') as file:
	for binid in sorted(ganon_map):
		print(ganon_map[binid], binid, file=file, sep="\t")

if overwrite:
	# move filter to ibf
	shutil.move(file_filter_in, file_ibf_out)
	# delete old files
	os.remove(file_nodes_in)
	os.remove(file_bins_in)
else:
	shutil.copy(file_filter_in, file_ibf_out)