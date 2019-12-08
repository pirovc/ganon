#!/usr/bin/env python3

import sys, pickle, gzip

db_prefix = sys.argv[1]

file_nodes = db_prefix+".nodes"
file_bins = db_prefix+".bins"
file_tax = db_prefix+".tax"
file_gnn = db_prefix+".gnn"

info = pickle.load(open(file_nodes, 'rb'))

tax = open(file_tax, "w")
print('\n'.join(['\t'.join([node,tax[0],tax[2],tax[1]]) for node,tax in info['filtered_nodes'].items()]), file=tax)
tax.close()

bins = []
with open(db_prefix+".bins", 'r') as file:
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
with gzip.open(file_gnn, 'wb') as f: pickle.dump(gnn, f)