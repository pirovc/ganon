#!/usr/bin/env python3
from collections import defaultdict

class Tax:
    def __init__(self, tax_files: list=None, ncbi_nodes: str=None, ncbi_names: str=None, ncbi_merged: str=None):
        # nodes[node] = (parent, rank, name)
        self.nodes = {}
        self.merged = {}
        # default root node
        self.nodes["1"] = ("0", "root", "root")
        if tax_files is not None:
            self.parse(tax_files)
        elif ncbi_nodes is not None:
            self.parse_ncbi_nodes(ncbi_nodes, ncbi_names)
            if ncbi_merged: self.parse_ncbi_merged(ncbi_merged)

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

    def parse_ncbi_nodes(self, ncbi_nodes, ncbi_names):
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

    def parse_ncbi_merged(self, ncbi_merged):
        self.merged = {}
        with open(ncbi_merged,'r') as file:
            for line in file:
                old_taxid, _, new_taxid, _ = line.split('\t',3)
                self.merged[old_taxid] = new_taxid

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

    def get_node(self, node):
        if node in self.nodes:
            return {"parent": self.nodes[node][0], "rank": self.nodes[node][1], "name": self.nodes[node][2]} 
        else:
            if node in self.merged and self.merged[node] in self.nodes:
                m = self.merged[node]
                return {"parent": self.nodes[m][0], "rank": self.nodes[m][1], "name": self.nodes[m][2]} 
            else:
                # If node not found, return root as parent
                return {"parent": "1", "rank": "na", "name": "na"}

    def get_node_rank_fixed(self, node, fixed_ranks):
        n = self.get_node(node)
        rank = n["rank"]
        while node!="0":
            if rank in fixed_ranks:
                return node, rank
            node = n["parent"]
            n = self.get_node(node)
            rank = n["rank"]
        return node, "na" #no standard rank identified
