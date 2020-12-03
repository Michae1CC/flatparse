import os

from flatparse.sequencing import FlatFileCreator

gff_path = os.path.join(os.getcwd(), 'test',
                        'data', 'Breviolum_minutum_short.gff')
output_path = os.path.join(os.getcwd(), 'test',
                           'data', 'tmp', 'test_flatfile_short.txt')

ffc_test = FlatFileCreator("SYMB1", "teamcx", gff_path, output_path)
ffc_test.create_flatfile()
