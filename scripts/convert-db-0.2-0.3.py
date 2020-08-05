#!/usr/bin/env python3
import os, sys, pickle, gzip, shutil

if len(sys.argv)==1:
	print("Conversion script for ganon database files from version 0.2.X to 0.3.X")
	print("Usage: ./convert-db-0.2-0.3.py db-prefix-in [db-prefix-out]")
	print("If no db-prefix-out is provided, the input database will be overwritten")
	sys.exit()

db_prefix_in = sys.argv[1]
overwrite=False
if len(sys.argv)==3:
	db_prefix_out = sys.argv[2]
else:
	db_prefix_out = db_prefix_in
	overwrite=True

file_tax_in = db_prefix_in+".tax"
file_gnn_in = db_prefix_in+".gnn"
file_map_in = db_prefix_in+".map"
file_ibf_in = db_prefix_in+".ibf"
file_tax_out = db_prefix_out+".tax"
file_gnn_out = db_prefix_out+".gnn"
file_map_out = db_prefix_out+".map"
file_ibf_out = db_prefix_out+".ibf"

gnn = pickle.load(gzip.open(file_gnn_in, 'rb'))

new_bins = []
for entry in gnn["bins"]:
	fields = entry.split("\t")

	# assembly field always at the end or empty
	if gnn["rank"]=="assembly":
		assembly = fields.pop(3)
		fields.append(assembly)
	else:
		fields.append("")

	# store fragment always, not with id
	if gnn["fragment_length"]>0:
		uid = fields.pop(0)
		i, pos = uid.split("/")
		st, en = pos.split(":")
		fields.insert(0, i)
		fields.insert(1, st)
		fields.insert(2, en)
	else:
		l = fields[1]
		fields.insert(1, "1")
		fields.insert(2, l)

	new_bins.append(fields)

gnn["bins"] = new_bins
with gzip.open(file_gnn_out, 'wb') as f: pickle.dump(gnn, f)

if not overwrite:
	shutil.copy(file_ibf_in, file_ibf_out)
	shutil.copy(file_tax_in, file_tax_out)
	shutil.copy(file_map_in, file_map_out)
