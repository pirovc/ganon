#!/usr/bin/env python3
import pandas as pd

class SeqInfo:
    seq_info_colums=['seqid','length','taxid','specialization']
    # Preset order without specialation
    seqinfo = pd.DataFrame(columns=seq_info_colums[:-1])
            
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

    def parse_seq_info_file(self, seq_info_file, parse_specialization: bool=False):
        self.seqinfo = pd.read_csv(seq_info_file, sep='\t', header=None, skiprows=0, index_col=False, names=self.seq_info_colums if parse_specialization else self.seq_info_colums[:-1])
        self.drop_duplicates()

    def drop_duplicates(self):
        self.seqinfo.drop_duplicates('seqid', keep="first", inplace=True)

    def validate_specialization(self):
        # convert length col to numeric int 
        # this is not set before due to issues of merging dataframes with different datatypes
        self.seqinfo['length'] = self.seqinfo['length'].astype(int)
        if 'specialization' in self.seqinfo.columns:
            # checl for invalid specialization
            idx_null_spec = self.seqinfo.specialization.isnull()
            # each specialization can have one parent taxid
            # get unique tuples taxid specialization
            taxid_spec = self.seqinfo[['taxid', 'specialization']].drop_duplicates()
            # check for multiple
            idx_multi_parent_spec = self.seqinfo.specialization.isin(taxid_spec.specialization[taxid_spec.specialization.duplicated(keep=False)].unique())
            # merge indices
            idx_replace = idx_null_spec | idx_multi_parent_spec
            if idx_replace.any():
                # replace duplicates with own seqid
                self.seqinfo.loc[idx_replace,"specialization"] = self.seqinfo.loc[idx_replace,"seqid"]
                return sum(idx_replace)
        return 0
