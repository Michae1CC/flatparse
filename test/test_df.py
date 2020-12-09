#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import os
import pandas as pd

example_table_path = os.path.join(os.getcwd(
), 'test', 'data', 'Slin_CCMP2456', 'S.linucheae_CCMP2456_uniprot_annotated.tsv')

example_table = pd.read_csv(
    example_table_path, engine='c', sep='\t', index_col='sequence')

print(example_table.head())
# print(example_table.loc['Slin_CCMP2456.gene6648.mRNA1']['accession'])
