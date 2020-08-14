import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import sys

# Usage plot_filter_size_dist.py acc_len_taxid.txt [nodes.dmp]
points=100
overlap_len=300
max_fp=0.05
hash_functions=3
kmer_size=19
max_size=0 #0 for default (1.5x min), in MB
rank="species"

def parse_input(file):  
    colums=['seqid', 'length', 'taxid', 'specialization']
    types={'seqid': 'str', 'length': 'uint64', 'taxid': 'str', 'specialization': 'str'}
    seqinfo = pd.read_table(file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)
    return seqinfo

def ibf_size_mb(bin_len, n_bins):
    return (math.ceil(-(1/((1-max_fp**(1/float(hash_functions)))**(1/float(hash_functions*(bin_len-kmer_size+1)))-1)))*optimal_bins(n_bins))/8388608

def approx_n_bins(bin_len): 
    frag_len=bin_len-overlap_len
    n_bins = sum([math.ceil(math.ceil(l/(frag_len-overlap_len))/(bin_len/(frag_len+overlap_len))) for l in groups_len.values()])
    return n_bins

def optimal_bins(n): 
    #return optimal number of bins for the IBF (multiples of 64)
    return (math.floor(n/64)+1)*64 

def parse_nodes_ranks(nodes_file):
    # READ nodes -> fields (1:TAXID 2:PARENT_TAXID 3:RANK)
    nodes = {}
    ranks = {}
    with open(nodes_file,'r') as fnodes:
        for line in fnodes:
            taxid, parent_taxid, rank, _ = line.split('\t|\t',3)
            taxid = taxid
            ranks[taxid] = rank
            nodes[taxid] = parent_taxid
    return nodes, ranks

def get_rank(ranks, node, rank):
    t = node
    try:
        while ranks[t]!=rank and t!="1": t = ranks[t]
    except:
        return node
    return t if t!="1" else node


seqinfo=parse_input(sys.argv[1])

if rank=="assembly":
    groups_len = seqinfo.groupby('assembly').sum().to_dict()['length']
elif rank=="taxid":
    groups_len = seqinfo.groupby('taxid').sum().to_dict()['length']
elif len(sys.argv)==3:
    _, ranks = parse_nodes_ranks(sys.argv[2])
    groups_len = pd.concat([seqinfo['taxid'].apply(lambda x: get_rank(ranks, x, rank)), seqinfo['length']], axis=1).groupby('taxid').sum().to_dict()['length']
else:
    print("Error parsing input")
    sys.exit(0)

max_bin_len = max(groups_len.values())
min_bin_len = 500

# %% Define time spans, initial values, and constants
bin_lens = np.geomspace(min_bin_len, max_bin_len, num=points)
n_bins = np.array([approx_n_bins(b) for b in bin_lens])
filter_sizes = np.array([ibf_size_mb(b, n_bins[i]) for i,b in enumerate(bin_lens)])

# keep only positives
pos_idx = filter_sizes>0
bin_lens = bin_lens[pos_idx]
n_bins = n_bins[pos_idx]
filter_sizes = filter_sizes[pos_idx]

min_size=filter_sizes.min()
if max_size==0: max_size = min_size*1.5
print("min. size: " + str(min_size) + "MB")
print("max. size: " + str(max_size) + "MB")
# keep only below max
idx_below_max = filter_sizes<=max_size

# If more than one valid point
if sum(idx_below_max)>1: 
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(bin_lens, filter_sizes, 'k--s')
    ax1.axhline(y=max_size, color='g', linestyle='--')
    ax1.axvline(x=bin_lens[idx_below_max][-1], color='g', linestyle='--')
    ax1.set_xlabel('Bin Length (bp)')
    ax1.set_ylabel('App. Filter Size (MB)')
    ax1b = ax1.twinx()
    ax1b.plot(bin_lens, n_bins , 'k--s', color="red")
    ax1b.set_ylabel('App. Number of bins', color="red")

    bin_lens = bin_lens[idx_below_max]
    n_bins = n_bins[idx_below_max]
    filter_sizes = filter_sizes[idx_below_max]
    # points between min and max to define good fit
    bin_lens = np.linspace(bin_lens.min(), bin_lens.max(), num=points)
    n_bins = np.array([approx_n_bins(b) for b in bin_lens])
    filter_sizes = np.array([ibf_size_mb(b, n_bins[i]) for i,b in enumerate(bin_lens)])
else:
    print("max size too small")
    sys.exit(0)

print("95 perc:")
print("bin_len:" + str(np.percentile(bin_lens, 95))+"bp")
print("size:" + str(np.percentile(filter_sizes, 95))+"MB")
print("n_bins:" + str(np.percentile(n_bins, 95)))

print("lowest n_bins:")
idx_min = np.where(n_bins == np.amin(n_bins))[0][0]
print(idx_min)
print("bin_len:" + str(bin_lens[idx_min]) +"bp")
print("size:" + str(filter_sizes[idx_min]) +"MB")
print("n_bins:" + str(n_bins[idx_min]))

ax2.plot(bin_lens, filter_sizes, 'k--s')
ax2.set_xlabel('Bin Length (bp)')
ax2.set_ylabel('App. Filter Size (MB)')
ax2b = ax2.twinx()
ax2b.plot(bin_lens, n_bins , 'k--s', color="red")
ax2b.set_ylabel('App. Number of bins', color="red")
plt.show()