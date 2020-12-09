#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import argparse

from flatparse.sequencing import FlatFileCreator

"""
Example usage:
    python -m flatparse.flatparse --locus_prefix="SYMB1" --dname="teamcx" --gff_path="D:\2020_SS\BioInfo\flatparse\test\data\test_data\Breviolum_minutum_short_v2.gff" --anno_path="D:\2020_SS\BioInfo\flatparse\test\data\test_data\Bminutum_top_uniprot_hits_v2.tsv" --output_path="D:\2020_SS\BioInfo\flatparse\test\data\tmp\test_flat.txt" --anno_delim=\t
"""


def run_ffc(args):

    main_ffc = FlatFileCreator(
        args.locus_prefix, args.dname, args.gff_path, args.anno_path,
        output_path=args.output_path, anno_delim=args.anno_delim)
    main_ffc.create_flatfile()

    return


def main():

    parser = argparse.ArgumentParser(description="Creates a flat file from a "
                                     "given gff file and annotations file.")

    parser.add_argument('--locus_prefix', type=str, required=True,
                        help='The prefix of each locus. Example: SYMB1 would be the prefix for the locus SYMB1_0001.')
    parser.add_argument('--dname', type=str, required=True,
                        help='The teamname that processed the data.')
    parser.add_argument('--gff_path', type=str, required=True,
                        help='A file path to the gff file.')
    parser.add_argument('--anno_path', type=str, required=True,
                        help="A file path to the annotation file. "
                        "Every annotation MUST have a column named 'sequence' "
                        "which is a unique identifier. Other optional columns "
                        "to include are: 'accession' 'function' 'gene_name' "
                        "'gene_synonym' and 'EC_number'. NOTE: these column "
                        "names are case sensitive.")

    parser.add_argument('--output_path', type=str, required=False, default=None,
                        help='A file path to output the contents of the flatfile. '
                        'Default output file is stdout.')
    parser.add_argument('--anno_delim', type=str,
                        required=False, default='\t',
                        help='A delimiter value for the annotation file.')

    args = parser.parse_args()
    run_ffc(args)

    exit(0)


if __name__ == "__main__":
    main()
