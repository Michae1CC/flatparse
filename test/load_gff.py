import os
from pprint import pprint
from flatparse.gff3 import Gff3

gff_path = gff_path = os.path.join(os.getcwd(), 'test',
                                   'data', 'test_data', 'Breviolum_minutum_short_v2.gff')
gff_path = gff_path = os.path.join(os.getcwd(), 'test',
                                   'data', 'Slin_CCMP2456', 'S.linucheae_CCMP2456_eg1.gff')

# initialize a Gff3 object
gff = Gff3(gff_file=gff_path)

pprint(len(gff.lines))
# pprint(gff.lines[2].keys())
pprint(gff.lines[0]['attributes']["Name"] + ".mRNA1")
# pprint(gff.lines[2]['attributes']['Parent'][0])
# pprint(gff.lines[1]['attributes'])
# pprint(gff.lines[1]['line_index'])
# pprint(gff.lines[1]['line_raw'])
# pprint(gff.lines[1]['start'])
# pprint(gff.lines[1]['end'])
# pprint(gff.lines[1]['type'])
