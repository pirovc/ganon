#!/usr/bin/env python3
import os, sys, pickle, gzip, shutil

if len(sys.argv)==1:
	print("Conversion script for ganon database files (only .gnn will change) from version 0.3.X to 0.4.X")
	print("Usage: ./convert-db-0.3-0.4.py db-prefix-in [db-prefix-out]")
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

#Add specialization and version fields
gnn["specialization"] = ""
gnn["version"] = "0.4.0"

# New nomenclature 
if gnn["rank"] == "taxid":
	gnn["rank"] = "leaves"
elif gnn["rank"] == "assembly":
	gnn["rank"] = "leaves"
	gnn["specialization"] = "assembly"

# write new .gnn
with gzip.open(file_gnn_out, 'wb') as f: pickle.dump(gnn, f)

if not overwrite:
	shutil.copy(file_ibf_in, file_ibf_out)
	shutil.copy(file_tax_in, file_tax_out)
	shutil.copy(file_map_in, file_map_out)

print("Conversion finished successfuly. Files written with the prefix: " + db_prefix_out)