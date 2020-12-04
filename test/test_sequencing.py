#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import os
import unittest

import flatparse.errors as errors
import flatparse.sequencing as sequencing
from flatparse.gff3 import Gff3

"""
Run:
    python -m unittest test.test_sequencing
with the venv activated to run all the tests within the test_sequencing script.
"""


class TestProteinInformation(unittest.TestCase):

    def test_completed_basic(self):

        simple_protein_info = sequencing.ProteinInformation(
            product="PROTEIN_ANNOTATION",
            teamname="teamcx",
            gene_locus="SYMB1_0001",
            sequence_prefix="mrna"
        )

        expected = "product PROTEIN_ANNOTATION\nprotein_id gnl|teamcx|SYMB1_0001\ntranscript_id gnl|teamcx|mrna.SYMB1_0001"

        self.assertEqual(str(simple_protein_info), expected)

    def test_incomplete_basic(self):

        simple_protein_info = sequencing.ProteinInformation(
            product=None,
            teamname="teamcx",
            gene_locus="SYMB1_0001",
            sequence_prefix="mrna"
        )

        try:
            str(simple_protein_info)
        except errors.InsufficientProteinInfo:
            # Correct error
            pass
        except Exception:
            self.fail()
        else:
            self.fail()

        simple_protein_info = sequencing.ProteinInformation(
            product="PROTEIN_ANNOTATION",
            teamname=None,
            gene_locus="SYMB1_0001",
            sequence_prefix="mrna"
        )

        try:
            str(simple_protein_info)
        except errors.InsufficientProteinInfo:
            # Correct error
            pass
        except Exception:
            self.fail()
        else:
            self.fail()

        simple_protein_info = sequencing.ProteinInformation(
            product="PROTEIN_ANNOTATION",
            teamname="teamcx",
            gene_locus=None,
            sequence_prefix="mrna"
        )

        try:
            str(simple_protein_info)
        except errors.InsufficientProteinInfo:
            # Correct error
            pass
        except Exception:
            self.fail()
        else:
            self.fail()

        simple_protein_info = sequencing.ProteinInformation(
            product="PROTEIN_ANNOTATION",
            teamname="teamcx",
            gene_locus="SYMB1_0001",
            sequence_prefix=None
        )

        try:
            str(simple_protein_info)
        except errors.InsufficientProteinInfo:
            # Correct error
            pass
        except Exception:
            self.fail()
        else:
            self.fail()

        simple_protein_info = sequencing.ProteinInformation(
            product=None,
            teamname=None,
            gene_locus=None,
            sequence_prefix=None
        )

        try:
            str(simple_protein_info)
        except errors.InsufficientProteinInfo:
            # Correct error
            pass
        except Exception:
            self.fail()
        else:
            self.fail()


class TestCDSGene(unittest.TestCase):

    # A file only containing two genes
    short_gff_path = os.path.join(os.getcwd(), 'test',
                                  'data', 'Breviolum_minutum_short.gff')
    short_gff = Gff3(short_gff_path)

    # A file only containing one gene
    gff_path_1 = os.path.join(os.getcwd(), 'test',
                              'data', 'test_data', 'Breviolum_minutum_1_gene.gff')
    gff_1 = Gff3(gff_path_1)

    loc_list_short_1 = [
        (26430, 26521),
        (26063,	26194),
        (24954,	25024),
        (24762,	24793),
        (24249,	24317),
        (23959,	24018),
        (23509,	23561),
        (22774,	22831),
        (22328,	22372),
        (21822,	21898),
        (21483,	21519),
        (20937,	21011),
        (20573,	20619),
        (20218,	20286),
        (19490,	19561),
        (19104,	19131),
        (18541,	18612),
        (18229,	18281),
        (9972,	10029),
        (9178,	9246),
        (8612,	8643),
        (8189,	8248),
        (7664,	7720),
        (7448,	7471),
        (6431,	6498),
        (6089,	6136),
        (5616,	5668)
    ]

    loc_list_short_2 = [
        (81761, 81789),
        (81221, 81289),
        (80642, 80733),
        (80096, 80213),
        (79546, 79686),
        (78969, 79092),
        (78371, 78478),
        (78076, 78103),
        (77459, 77546),
        (77115, 77161),
        (76913, 76962),
        (76678, 76803),
        (75787, 75858),
        (75505, 75581),
        (75263, 75341),
        (75024, 75127),
        (74704, 74752),
        (74174, 74269),
        (73783, 73943),
        (73228, 73317),
        (72776, 72861),
        (71780, 71871),
        (71341, 71470),
        (70762, 70870),
        (70016, 70069),
        (69728, 69783),
        (69386, 69450),
        (67007, 67114),
        (66502, 66549),
        (65882, 65946),
        (65171, 65256),
        (64603, 64669),
        (64142, 64224),
        (63695, 63763),
        (63463, 63488),
        (62882, 62978),
        (62494, 62546),
        (61815, 61878),
        (61605, 61640),
        (61186, 61326),
        (60457, 60553),
        (60108, 60177),
        (59566, 59628),
        (59213, 59269),
        (58818, 58853),
        (58160, 58200),
        (57578, 57663),
        (56681, 56756),
        (56137, 56220),
        (55693, 55805),
        (55140, 55182),
        (54696, 54792),
        (54404, 54447),
        (53720, 53825)
    ]

    def test_construction(self):

        cls = self.__class__

        sequencing.CDSGene("SYMB1", "teamcx", cls.short_gff, 1, 3, 55)

    def test_simple_get_lbp_list(self):

        cls = self.__class__

        cds_simple = sequencing.CDSGene(
            "SYMB1", "teamcx", cls.short_gff, 1, 3, 55)

        self.assertEqual(cds_simple.get_lbp_list(), cls.loc_list_short_1)

        cds_simple = sequencing.CDSGene(
            "SYMB1", "teamcx", cls.short_gff, 2, 59, 165)

        self.assertEqual(cds_simple.get_lbp_list(), cls.loc_list_short_2)

    def test_1_gene_get_lbp_list(self):

        cls = self.__class__

        cds_simple = sequencing.CDSGene(
            "SYMB1", "teamcx", cls.gff_1, 1, 3, 55)

        self.assertEqual(cds_simple.get_lbp_list(), cls.loc_list_short_1)

    def test_simple_str(self):

        cls = self.__class__

        cds_simple = sequencing.CDSGene(
            "SYMB1", "teamcx", cls.short_gff, 1, 3, 55)

        first_line = "26430\t26521\tCDS"

        intermediary_locs = '\n'.join(
            map(lambda tup_: "{}\t{}".format(tup_[0], tup_[1]), cls.loc_list_short_1[1:]))

        protein_info = "\t\t\tproduct PROTEIN_ANNOTATION\n" + \
            "\t\t\tprotein_id gnl|teamcx|SYMB1_0001\n" + \
            "\t\t\ttranscript_id gnl|teamcx|mrna.SYMB1_0001"

        expected = first_line + '\n' + intermediary_locs + '\n' + protein_info

        self.assertEqual(str(cds_simple), expected)

        # Test the second gene

        cds_simple = sequencing.CDSGene(
            "SYMB1", "teamcx", cls.short_gff, 2, 59, 165)

        first_line = "81761\t81789\tCDS"

        intermediary_locs = '\n'.join(
            map(lambda tup_: "{}\t{}".format(tup_[0], tup_[1]), cls.loc_list_short_2[1:]))

        protein_info = "\t\t\tproduct PROTEIN_ANNOTATION\n" + \
            "\t\t\tprotein_id gnl|teamcx|SYMB1_0002\n" + \
            "\t\t\ttranscript_id gnl|teamcx|mrna.SYMB1_0002"

        expected = first_line + '\n' + intermediary_locs + '\n' + protein_info

        self.assertEqual(str(cds_simple), expected)


class TestmRNAGene(unittest.TestCase):

    short_gff = Gff3(os.path.join(os.getcwd(), 'test',
                                  'data', 'Breviolum_minutum_short.gff'))

    # A file only containing one gene
    gff_path_1 = os.path.join(os.getcwd(), 'test',
                              'data', 'test_data', 'Breviolum_minutum_1_gene.gff')
    gff_1 = Gff3(gff_path_1)

    loc_list_short_1 = [
        (5616, 26521),
    ]

    loc_list_short_2 = [
        (53720, 81789),
    ]

    def test_construction(self):

        cls = self.__class__

        sequencing.mRNAGene("SYMB1", "teamcx", cls.short_gff, 1, 1, 2)

    def test_simple_get_lbp_list(self):

        cls = self.__class__

        mRNA_simple = sequencing.mRNAGene(
            "SYMB1", "teamcx", cls.short_gff, 1, 1, 2)

        self.assertEqual(mRNA_simple.get_lbp_list(), cls.loc_list_short_1)

        mRNA_simple = sequencing.mRNAGene(
            "SYMB1", "teamcx", cls.short_gff, 2, 57, 58)

        self.assertEqual(mRNA_simple.get_lbp_list(), cls.loc_list_short_2)

        # Check to see if single line mRNA won't break
        mRNA_simple = sequencing.mRNAGene(
            "SYMB1", "teamcx", cls.short_gff, 2, 57, 57)

        self.assertEqual(mRNA_simple.get_lbp_list(), cls.loc_list_short_2)

    def test_1_gene_get_lbp_list(self):

        cls = self.__class__

        mRNA_simple = sequencing.mRNAGene(
            "SYMB1", "teamcx", cls.gff_1, 1, 1, 2)

        self.assertEqual(mRNA_simple.get_lbp_list(), cls.loc_list_short_1)

        mRNA_simple = sequencing.mRNAGene(
            "SYMB1", "teamcx", cls.short_gff, 2, 57, 58)

        self.assertEqual(mRNA_simple.get_lbp_list(), cls.loc_list_short_2)

    def test_simple_str(self):

        cls = self.__class__

        mRNA_simple = sequencing.mRNAGene(
            "SYMB1", "teamcx", cls.short_gff, 1, 1, 2)

        first_line = "5616\t26521\tmRNA"

        protein_info = "\t\t\tproduct PROTEIN_ANNOTATION\n" + \
            "\t\t\tprotein_id gnl|teamcx|SYMB1_0001\n" + \
            "\t\t\ttranscript_id gnl|teamcx|mrna.SYMB1_0001"

        expected = first_line + '\n' + protein_info

        self.assertEqual(str(mRNA_simple), expected)

        mRNA_simple = sequencing.mRNAGene(
            "SYMB1", "teamcx", cls.short_gff, 2, 57, 58)

        first_line = "53720\t81789\tmRNA"

        protein_info = "\t\t\tproduct PROTEIN_ANNOTATION\n" + \
            "\t\t\tprotein_id gnl|teamcx|SYMB1_0002\n" + \
            "\t\t\ttranscript_id gnl|teamcx|mrna.SYMB1_0002"

        expected = first_line + '\n' + protein_info

        self.assertEqual(str(mRNA_simple), expected)

    def test_1_gene_str(self):

        cls = self.__class__

        mRNA_simple = sequencing.mRNAGene(
            "SYMB1", "teamcx", cls.gff_1, 1, 1, 2)

        first_line = "5616\t26521\tmRNA"

        protein_info = "\t\t\tproduct PROTEIN_ANNOTATION\n" + \
            "\t\t\tprotein_id gnl|teamcx|SYMB1_0001\n" + \
            "\t\t\ttranscript_id gnl|teamcx|mrna.SYMB1_0001"

        expected = first_line + '\n' + protein_info

        self.assertEqual(str(mRNA_simple), expected)


class TestFeature(unittest.TestCase):

    # A file only containing two genes
    short_gff_path = os.path.join(os.getcwd(), 'test',
                                  'data', 'Breviolum_minutum_short.gff')
    short_gff = Gff3(short_gff_path)

    # A file only containing one gene
    gff_path_1 = os.path.join(os.getcwd(), 'test',
                              'data', 'test_data', 'Breviolum_minutum_1_gene.gff')
    gff_1 = Gff3(gff_path_1)

    def test_construction(self):

        cls = self.__class__

        sequencing.GeneSequence(
            "SYMB1", "teamcx", cls.short_gff, 1, "Bmin.gene1", 0, 56
        )

    def test_simple_analyse_sequence(self):

        cls = self.__class__

        test_gene_seq = sequencing.GeneSequence(
            "SYMB1", "teamcx", cls.short_gff, 1, "Bmin.gene1", 0, 55
        )

        expected_sequence = {
            (1, 2): "mRNA",
            (3, 55): "CDS"
        }

        self.assertEqual(test_gene_seq.analyse_sequence(), expected_sequence)

        # Test the second gene
        test_gene_seq = sequencing.GeneSequence(
            "SYMB1", "teamcx", cls.short_gff, 2, "Bmin.gene1", 56, 165
        )

        expected_sequence = {
            (57, 58): "mRNA",
            (59, 165): "CDS"
        }

        self.assertEqual(test_gene_seq.analyse_sequence(), expected_sequence)

    def test_1_gene_analyse_sequence(self):

        cls = self.__class__

        test_gene_seq = sequencing.GeneSequence(
            "SYMB1", "teamcx", cls.gff_1, 1, "Bmin.gene1", 0, 55
        )

        expected_sequence = {
            (1, 2): "mRNA",
            (3, 55): "CDS"
        }

        self.assertEqual(test_gene_seq.analyse_sequence(), expected_sequence)

    def test_simple_str(self):

        cls = self.__class__

        test_gene_seq = sequencing.GeneSequence(
            "SYMB1", "teamcx", cls.short_gff, 1, "Bmin.gene1", 0, 56
        )

        mRNA_simple = sequencing.mRNAGene(
            "SYMB1", "teamcx", cls.short_gff, 1, 1, 2)

        cds_simple = sequencing.CDSGene(
            "SYMB1", "teamcx", cls.short_gff, 1, 3, 55)

        line_list = [
            ">Features Bmin.gene1",
            "5616\t26521\tgene",
            "\t\t\tlocus_tag\tSYMB1_0001",
        ]

        line_list.append(str(mRNA_simple))
        line_list.append(str(cds_simple))

        expected = '\n'.join(line_list)

        self.assertEqual(str(test_gene_seq), expected)

        # Test the second gene
        test_gene_seq = sequencing.GeneSequence(
            "SYMB1", "teamcx", cls.short_gff, 2, "Bmin.gene2", 56, 165
        )

        mRNA_simple = sequencing.mRNAGene(
            "SYMB1", "teamcx", cls.short_gff, 2, 57, 58)

        cds_simple = sequencing.CDSGene(
            "SYMB1", "teamcx", cls.short_gff, 2, 59, 165)

        line_list = [
            ">Features Bmin.gene2",
            "53720\t81789\tgene",
            "\t\t\tlocus_tag\tSYMB1_0002",
        ]

        line_list.append(str(mRNA_simple))
        line_list.append(str(cds_simple))

        expected = '\n'.join(line_list)

        self.assertEqual(str(test_gene_seq), expected)

    def test_1_gene_str(self):

        cls = self.__class__

        test_gene_seq = sequencing.GeneSequence(
            "SYMB1", "teamcx", cls.gff_1, 1, "Bmin.gene1", 0, 55
        )

        mRNA_simple = sequencing.mRNAGene(
            "SYMB1", "teamcx", cls.short_gff, 1, 1, 2)

        cds_simple = sequencing.CDSGene(
            "SYMB1", "teamcx", cls.short_gff, 1, 3, 55)

        line_list = [
            ">Features Bmin.gene1",
            "5616\t26521\tgene",
            "\t\t\tlocus_tag\tSYMB1_0001",
        ]

        line_list.append(str(mRNA_simple))
        line_list.append(str(cds_simple))

        expected = '\n'.join(line_list)

        self.assertEqual(str(test_gene_seq), expected)


class TestFlatFileCreator(unittest.TestCase):

    # A file only containing two genes
    short_gff_path = os.path.join(os.getcwd(), 'test',
                                  'data', 'Breviolum_minutum_short.gff')
    short_gff = Gff3(short_gff_path)

    # A file only containing one gene
    gff_path_1 = os.path.join(os.getcwd(), 'test',
                              'data', 'test_data', 'Breviolum_minutum_1_gene.gff')
    gff_1 = Gff3(gff_path_1)

    def test_constructor(self):

        cls = self.__class__

        sequencing.FlatFileCreator("SYMB1", "teamcx", cls.short_gff_path)

    def test_constructor_invalid_gff_path(self):

        invalid_path_ext = os.path.join(os.getcwd(), 'test',
                                        'data', 'test_data', 'Breviolum_minutum.txt')

        non_existant_path = os.path.join(os.getcwd(), 'test',
                                         'data', 'test_data', 'Breviolum_minutum.gff')

        try:
            sequencing.FlatFileCreator("SYMB1", "teamcx", invalid_path_ext)
        except OSError:
            # Correct error
            pass
        except Exception:
            self.fail()
        else:
            self.fail()

        try:
            sequencing.FlatFileCreator("SYMB1", "teamcx", non_existant_path)
        except OSError:
            # Correct error
            pass
        except Exception:
            self.fail()
        else:
            self.fail()

    def test_simple_get_gene_dict(self):

        cls = self.__class__

        test_ffc = sequencing.FlatFileCreator(
            "SYMB1", "teamcx", cls.short_gff_path)

        # Expected output
        expected = {
            "Bmin.gene1": (0, 55),
            "Bmin.gene2": (56, 165)
        }

        self.assertEqual(test_ffc.get_gene_dict(), expected)

    def test_1_gene_get_gene_dict(self):

        cls = self.__class__

        test_ffc = sequencing.FlatFileCreator(
            "SYMB1", "teamcx", cls.short_gff_path)

        # Expected output
        expected = {
            "Bmin.gene1": (0, 55),
            'Bmin.gene2': (56, 165)
        }

        self.assertEqual(test_ffc.get_gene_dict(), expected)

    def test_single_get_gene_dict(self):

        cls = self.__class__

        test_ffc = sequencing.FlatFileCreator(
            "SYMB1", "teamcx", cls.gff_path_1)

        # Expected output
        expected = {
            "Bmin.gene1": (0, 55)
        }

        self.assertEqual(test_ffc.get_gene_dict(), expected)

    def test_simple_str(self):

        cls = self.__class__

        gene_seq_1 = sequencing.GeneSequence(
            "SYMB1", "teamcx", cls.short_gff, 1, "Bmin.gene1", 0, 55
        )

        gene_seq_2 = sequencing.GeneSequence(
            "SYMB1", "teamcx", cls.short_gff, 2, "Bmin.gene2", 56, 165
        )

        expected = str(gene_seq_1) + '\n' + str(gene_seq_2)

        ffc_test = sequencing.FlatFileCreator(
            "SYMB1", "teamcx", cls.short_gff_path)

        self.assertEqual(str(ffc_test), expected)

    def test_1_gene_str(self):

        cls = self.__class__

        gene_seq_1 = sequencing.GeneSequence(
            "SYMB1", "teamcx", cls.gff_1, 1, "Bmin.gene1", 0, 55
        )

        expected = str(gene_seq_1)

        ffc_test = sequencing.FlatFileCreator(
            "SYMB1", "teamcx", cls.gff_path_1)

        self.assertEqual(str(ffc_test), expected)


if __name__ == '__main__':
    unittest.main()
