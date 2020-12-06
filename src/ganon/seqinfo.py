#!/usr/bin/env python3
import pandas as pd
from io import StringIO
from ganon.util import *

class SeqInfo:
    seq_info_colums=['seqid','length','taxid','specialization']
    # pd.Int64Dtype() nullable INT https://pandas.pydata.org/pandas-docs/stable/user_guide/integer_na.html#integer-na
    seq_info_types={'seqid': 'str', 'length': 'int', 'taxid': 'str', 'specialization': 'str'}
    
    # Preset order without specialation
    seqinfo = pd.DataFrame(columns=seq_info_colums[:-1])

    #def __init__(self):
            
    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
        return 'SeqInfo({})'.format(', '.join(args))

    def to_csv(self):
        return self.seqinfo.to_csv(sep="\t", header=False, index=False)

    def size(self):
        return self.seqinfo.shape[0]

    def get_seqids(self):
        return self.seqinfo.seqid

    def remove_seqids(self, seqids):
        self.seqinfo = self.seqinfo[~self.seqinfo['seqid'].isin(seqids)]

    def write_seqid_file(self, seqid_file):
        self.seqinfo["seqid"].to_csv(seqid_file, sep="\t", header=False, index=False)

    def write(self, file):
        self.seqinfo.to_csv(file, sep="\t", header=False, index=False)

    def drop_zeros(self, col: str):
        self.seqinfo = self.seqinfo[self.seqinfo[col]>0]

    def append(self, df):
        self.seqinfo = self.seqinfo.append(df, ignore_index=True, sort=False)

    def clear(self):
        self.seqinfo = self.seqinfo[0:0]

    def paste_cols(self, cols, df):
        self.seqinfo[cols] = df        

    def dropna(self, subset):
        self.seqinfo.dropna(subset=subset, inplace=True)

    def join(self, df, field):
        self.seqinfo[field] = self.seqinfo.join(df.set_index('seqid'), on="seqid", how="left", rsuffix="_tojoin")[field+"_tojoin"]        

    def parse_seq_info_file(self, seq_info_file):
        self.seqinfo = pd.read_csv(seq_info_file, sep='\t', header=None, skiprows=0, names=self.seq_info_colums)

    def validate(self):
        self.seqinfo['length'] = self.seqinfo['length'].astype(int)