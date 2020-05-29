#!/usr/bin/env python3

import sys
import argparse
import pickle
import gzip
from io import StringIO
import pandas as pd
from collections import defaultdict

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument('-q', '--fastq-file', required=False, type=str, metavar='',  help='Input reads used on ganon classify (.fq[.gz])')
	parser.add_argument('-g', '--ganon-output', type=str, dest='ganon_output', required=True, help='Results from ganon classify (.lca or .all)')
	parser.add_argument('-d', '--db-prefix', required=True, type=str, dest='db_prefix', help='Database input prefix used on ganon classify')
	parser.add_argument('-o', '--output-prefix', type=str, dest='output_prefix', required=True, help='Output prefix')
	parser.add_argument('-i', '--input-files', required=False, type=str, nargs="*", metavar='',  help='Input reference sequence fasta files (.fna[.gz]) used on ganon build. If provided, creates a set of references for each output')
	parser.add_argument('-r', '--rank', type=str, default='species', dest="rank", help='Target taxonomic to output [assembly,taxid,species,genus,...]. It cannot be more specific than the rank used on ganon build. If assingment rank is less precise than rank chosen, assignment rank is reported. Default: species')
	args = parser.parse_args()

	bins = parse_bins(args.db_prefix)
	tax = parse_tax(args.db_prefix)
	# add row with target on the chosen rank
	bins["target"] = bins.apply(lambda row : get_rank(tax, row['group'], args.rank), axis = 1) 

	unique_targets = list(pd.unique(bins["target"]))
	
	# group is where the assignments are done, target is the choosen target output rank taxid
	group_target = dict(bins[["group","target"]].values)
	
	out = {}
	for target in unique_targets:
		out[target] = open(args.output_prefix + "." + tax[target][1] + "." + target + ".out", "w")
		 
	# Store all and write specific match
	am = defaultdict(set)
	with open(args.ganon_output) as file: 
		for line in file:
			field = line.rstrip().split("\t")
			target = group_target[field[1]]
			am[field[0]].add(target)
			print(field[0], field[1], file=out[target], sep="\t")
			
	# close out
	for target in unique_targets:
		out[target].close()
	

	fastq_out = {}
	for target in unique_targets:
		fastq_out[target] = open(args.output_prefix + "." + tax[target][1] + "." + target + ".fq", "w")
	for name, seq, qual in readfq(open(args.fastq_file, "r")):
		for target in am[name]:
			print("@"+name, seq, "+", qual, file=fastq_out[target], sep="\n")

	# close fastq
	for target in unique_targets:
		fastq_out[target].close()
	
	# if input references are provided
	if args.input_files:
		seqid_target = dict(bins[["seqid","target"]].values)
		
		#open fastq files
		fasta_out = {}
		for target in unique_targets:
			fasta_out[target] = open(args.output_prefix + "." + tax[target][1] + "." + target + ".fna", "w")
		
		for ref in args.input_files:
			file = open(ref, "r") if not ref.endswith(".gz") else gzip.open(ref, "rt")
			for name, seq, _ in readfq(file):
				if name in seqid_target:
					print(">"+name, seq, file=fasta_out[seqid_target[name]], sep="\n")
			file.close()

		# close fasta files
		for target in unique_targets:
			fasta_out[target].close()

	# print matching info
	for index, row in bins.iterrows():
		print(row["seqid"], row["group"], row["target"], sep="\t")

def parse_tax(db_prefix):
	tax = {}
	with open(db_prefix + ".tax" , 'r') as file:
		for line in file:
			node, parent, rank, name = line.rstrip().split("\t")
			if node not in tax:
				tax[node] = (parent,rank,name)
	return tax

def get_rank(tax, node, rank):
	t = node
	try:
		while tax[t][1]!=rank and t!="1": t = tax[t][0]
	except:
		return node
	return t if t!="1" else node

def parse_bins(db_prefix):
	gnn = pickle.load(gzip.open(db_prefix + ".gnn", "rb"))
	bins = pd.DataFrame([line.split("\t") for line in gnn["bins"]],
							columns=['seqid','length', 'taxid', 'binid'] if gnn["rank"]!="assembly" else ['seqid','length', 'taxid', 'group', 'binid'])
	if gnn["rank"]!="assembly":
		bins['group'] = bins['taxid'] # Add taxid as group
	if gnn["fragment_length"]: 
		bins[['seqid','begin','end']] = bins.seqid.str.split(r"\/|:", expand=True)
	else:
		bins['begin'] = 1
		bins['end'] = bins['length']
	return bins

def readfq(fp): # this is a generator function
	last = None # this is a buffer keeping the last unprocessed line
	while True: # mimic closure; is it a bad idea?
		if not last: # the first record or a record following a fastq
			for l in fp: # search for the start of the next record
				if l[0] in '>@': # fasta/q header line
					last = l[:-1] # save this line
					break
		if not last: break
		name, seqs, last = last[1:].partition(" ")[0], [], None
		for l in fp: # read the sequence
			if l[0] in '@+>':
				last = l[:-1]
				break
			seqs.append(l[:-1])
		if not last or last[0] != '+': # this is a fasta record
			yield name, ''.join(seqs), None # yield a fasta record
			if not last: break
		else: # this is a fastq record
			seq, leng, seqs = ''.join(seqs), 0, []
			for l in fp: # read the quality
				seqs.append(l[:-1])
				leng += len(l) - 1
				if leng >= len(seq): # have read enough quality
					last = None
					yield name, seq, ''.join(seqs); # yield a fastq record
					break
			if last: # reach EOF before reading enough quality
				yield name, seq, None # yield a fasta record instead
				break

if __name__ == '__main__':
	main()
