#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import os
import pandas as pd

anno_filepath = os.path.join(os.getcwd(), 'test',
                             'data', 'test_data', 'Bminutum_top_uniprot_hits_v2.tsv')
# anno_df = pd.read_csv(anno_filepath, sep='\t', header=None, index_col=0)
anno_df = pd.read_csv(anno_filepath, sep=None, header=None, index_col=0)

print(anno_df.head())
# Get all the annotation values belonging to the gene Bmin.gene1.mRNA1
print(anno_df.loc["Bmin.gene1.mRNA1"])
print(anno_df.loc["Bmin.gene1.mRNA1"][1])
