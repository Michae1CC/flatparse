import os

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


def main():
    ffc_cls()


if __name__ == '__main__':
    main()
