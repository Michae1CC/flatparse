import os
from pprint import pprint
from flatparse.gff3 import Gff3

gff_path = gff_path = os.path.join(os.getcwd(), 'test',
                                   'data', 'test_data', 'Breviolum_minutum_short_v2.gff')

# initialize a Gff3 object
gff = Gff3(gff_file=gff_path)

pprint(len(gff.lines))
pprint(gff.lines[0].keys())
pprint(gff.lines[0]['seqid'])
# pprint(gff.lines[1]['attributes'])
# pprint(gff.lines[1]['line_index'])
# pprint(gff.lines[1]['line_raw'])
# pprint(gff.lines[1]['start'])
# pprint(gff.lines[1]['end'])
# pprint(gff.lines[1]['type'])
