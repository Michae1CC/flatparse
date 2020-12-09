#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import re
import os
import pandas as pd

example_table_path = os.path.join(os.getcwd(
), 'test', 'data', 'Slin_CCMP2456', 'S.linucheae_CCMP2456_uniprot_annotated.tsv')

example_table = pd.read_csv(
    example_table_path, engine='c', sep='\t', index_col='sequence')

# print(example_table.head())
# print(example_table.loc['Slin_CCMP2456.gene6648.mRNA1']['accession'])
example_function = example_table.loc['Slin_CCMP2456.gene6648.mRNA1']['function']
example_function = example_table.loc['Slin_CCMP2456.gene18741.mRNA1']['function']
print(example_function)
# EC_pattern = r'\(EC (?P<ID>[0-9_.-]+)\)'
# # re.findall(EC_pattern, example_function)
# result = re.search(EC_pattern, example_function)
# EC_number = result.group("ID")
# print(EC_number)

# example_function = re.sub(EC_pattern, '', example_function)
# print(example_function)
