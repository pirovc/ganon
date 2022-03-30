#!/usr/bin/env python3
import pandas as pd

class Bins:
    # bins columns pandas dataframe
    columns=['seqid', 'seqstart', 'seqend', 'length', 'taxid', 'binid', 'specialization']

    def __init__(self, taxsbp_ret: list, use_specialization: bool=False):
        if taxsbp_ret:
            self.bins = pd.DataFrame(taxsbp_ret)
            # Always returns 7 cols from taxsbp, remove empty col if specialization not set
            if not use_specialization and self.bins.shape[1]==7: del self.bins[6]
            self.bins.rename(columns={i:n for i,n in enumerate(self.columns)}, inplace=True)
        else:
            self.bins = pd.DataFrame(columns=self.columns)

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
        return 'Bins({})'.format(', '.join(args))

    def size(self):
        return self.bins.shape[0]

    def get_subset(self,seqids: list=[], binids: list=[]):
        subBins = Bins()
        if seqids: subBins.bins = self.bins.loc[self.bins['seqid'].isin(seqids)]
        elif binids: subBins.bins = self.bins.loc[self.bins['binid'].isin(binids)]
        return subBins

    def to_csv(self):
        return self.bins.to_csv(sep="\t",header=False, index=False)

    def get_list(self):
        return self.bins.values.tolist()

    def get_bins(self):
        return self.bins

    def remove_seqids(self, seqids):
        self.bins = self.bins.loc[~self.bins['seqid'].isin(seqids)]

    def get_number_of_bins(self):
        return self.get_binids().size
    
    def get_taxids(self):
        return self.bins.taxid.unique()

    def get_specialization(self):
        return self.bins.specialization.unique()
        
    def get_binids(self):
        return self.bins.binid.unique()

    def get_seqids(self):
        return self.bins.seqid.unique()

    def get_binid_length_sum(self):
        return self.bins.groupby(['binid'])['length'].sum()
        
    def get_max_bin_length(self):
        return self.get_binid_length_sum().max()

    def get_max_split_bins(self, use_specialization: bool=False):
        if use_specialization:
            return self.bins[['specialization','binid']].drop_duplicates().value_counts('specialization').max()
        else:
            return self.bins[['taxid','binid']].drop_duplicates().value_counts('taxid').max()

    def get_specialization_taxid(self):
        return self.bins[['specialization','taxid']].drop_duplicates().set_index('specialization')

    def merge(self, bins):
        self.bins = pd.concat([self.bins, bins.get_bins()])

    def write_acc_bin_file(self, acc_bin_file, binids: set=None):
        if binids is None:
            self.bins.to_csv(acc_bin_file, header=False, index=False, columns=['seqid','seqstart','seqend','binid'], sep='\t')
        else:
            self.bins.loc[self.bins['binid'].isin(binids)].to_csv(acc_bin_file, header=False, index=False, columns=['seqid','seqstart','seqend','binid'], sep='\t')

    def write_map_file(self, map_file, use_specialization: bool=False):
        if use_specialization:
            self.bins[['specialization','binid']].drop_duplicates().to_csv(map_file,header=False, index=False, sep='\t')
        else:
            self.bins[['taxid','binid']].drop_duplicates().to_csv(map_file,header=False, index=False, sep='\t')