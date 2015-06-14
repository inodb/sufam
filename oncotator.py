"""
Functions for oncotator
"""
import pandas as pd
import pyhgvs as hgvs

class Oncotator(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.df = pd.read_csv(filepath, sep="\t", comment="#")

    def get_hgvs_mutations(self, transcript_id):
        return [hgvs.HGVSName(r) for r in self.df.HGVS_coding_DNA_change \
            if r != "Exception_encountered" and \
            hgvs.HGVSName(r).transcript == transcript_id]
