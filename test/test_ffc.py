import os

import argparse
import warnings
from pprint import pprint

from flatparse.gff3 import Gff3
from flatparse.sequencing import FlatFileCreator
from flatparse.sequencing import Feature

gff_path = os.path.join(os.getcwd(), 'test',
                        'data', 'test_data', 'Breviolum_minutum_short_v2.gff')
anno_path_v2 = os.path.join(os.getcwd(), 'test',
                            'data', 'test_data', 'Bminutum_top_uniprot_hits_v2.tsv')

output_path = os.path.join(os.getcwd(), 'test',
                           'data', 'tmp', 'test_flatfile_short.txt')

tmp_out_path = os.path.join(
    os.getcwd(), "test", "data", "tmp", "test_flat.txt")


def ffc_cls():

    ffc_test = FlatFileCreator(
        "SYMB1", "teamcx", gff_path, anno_path_v2, output_path=output_path)

    # pprint(ffc_test.analyse())

    ffc_test.create_flatfile()


def feature_cls():

    kwargs = {"dname": "teamcx", "locus_prefix": "SYMB1"}

    gff = Gff3(gff_file=gff_path)
    test_feat = Feature("Bmin.scaffold2", gff, 0, 901, **kwargs)

    # pprint(test_feat.analyse())
    print(str(test_feat)[:200])


def test_delim():

    parser = argparse.ArgumentParser(description="Creates a flat file from a "
                                     "given gff file and annotations file.")

    parser.add_argument('--delimiter', type=str,
                        required=False, default='\t')
    args = parser.parse_args()

    delimiter = args.delimiter

    for old, new in [('\\n', '\n'), ('\\t', '\t'), ('\\r', '\r')]:
        delimiter = delimiter.replace(old, new)

    warnings.warn('This is a warning only {0!r}'.format(
        delimiter), RuntimeWarning)

    print('Hello ' + delimiter + ' world')


def main():
    test_delim()


if __name__ == '__main__':
    main()
