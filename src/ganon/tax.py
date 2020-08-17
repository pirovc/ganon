#!/usr/bin/env python3
from collections import defaultdict

class Tax:
    def __init__(self, tax_files: list=None, ncbi_nodes: str=None, ncbi_names: str=None):
        # nodes[node] = (parent, rank, name)
        self.nodes = {}
        # default root node
        self.nodes["1"] = ("0", "no rank", "root")
        if tax_files is not None:
            self.parse(tax_files)
        elif ncbi_nodes is not None:
            self.parse_ncbi(ncbi_nodes, ncbi_names)

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
        return 'Tax({})'.format(', '.join(args))

    def add_nodes(self, new_nodes, label):
        # add nodes from a pandas dataframe [group taxid]
        for node,parent in new_nodes.iterrows():
            if node not in self.nodes:
                self.nodes[node] = (parent['taxid'],label,node) # repeat node on name

    def merge(self, extra_tax): # duplicates solved by current tax
        for node in extra_tax.nodes:
            if node not in self.nodes: # if present in the current tax
                self.nodes[node] = extra_tax.nodes[node]

    def filter(self, filter_nodes):
        new_nodes = {}
        # Filter nodes for used taxids
        for node in filter_nodes:
            if node in self.nodes: # if present in the current tax
                t = node
                while t!="0": # while not at root, get lineage
                    if t in new_nodes: break # branch already added
                    new_nodes[t] = self.nodes[t] # copy node
                    t = self.nodes[t][0] # get parent
        self.nodes = new_nodes

    def parse(self, tax_files):
        for tax_file in tax_files:
            with open(tax_file , 'r') as file:
                for line in file:
                    node, parent, rank, name = line.rstrip().split("\t")
                    if node not in self.nodes:
                        self.nodes[node] = (parent,rank,name)

    def parse_ncbi(self, ncbi_nodes, ncbi_names):
        names = defaultdict(lambda: "")
        if ncbi_names is not None:
            with open(ncbi_names,'r') as file:
                for line in file:
                    node, name, _, name_class = line.split('\t|\t') # READ names -> fields (1:TAXID 2:NAME 3:UNIQUE NAME 4:NAME CLASS)
                    if name_class.replace('\t|\n','')=="scientific name":
                        names[node] = name
        with open(ncbi_nodes , 'r') as file:
            for line in file:
                node, parent, rank, _ = line.split('\t|\t', 3) # READ nodes -> fields (1:TAXID 2:PARENT_TAXID 3:RANK)
                if node not in self.nodes:
                    self.nodes[node] = (parent,rank,names[node])

    def write(self, file):
        # .tax: taxid/assembly <tab> parent taxid <tab> rank <tab> name
        tax_file = open(file,'w')
        for node,(parent, rank, name) in self.nodes.items():
            print(node, parent, rank, name, sep="\t", file=tax_file)
        tax_file.close()

    def get_rank(self, node, rank):
        t = node
        try:
            while self.nodes[t][1]!=rank and t!="1": t = self.nodes[t][0]
        except:
            return node
        return t if t!="1" else node

    def get_node_rank_fixed(self, node, fixed_ranks):
        if node!="0":
            original_rank = self.nodes[node][1]
            original_taxid = node
            while node!="0":
                if(self.nodes[node][1] in fixed_ranks):
                    #everything below species (not being assembly) is counted as species+
                    if "species+" in fixed_ranks and original_rank!="species" and original_rank!="assembly" and self.nodes[node][1]=="species":
                        return original_taxid, "species+"
                    else:
                        return node, self.nodes[node][1]
                node = self.nodes[node][0]
            return "1", "root" #no standard rank identified
        else:
            return "0", ""
