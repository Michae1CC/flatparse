#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import abc
import warnings
import pandas as pd
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple, Type, Union

import flatparse.errors as errors
import flatparse.path as path
from flatparse.gff3 import Gff3


class ProteinInformation:
    """
    Holds information about the protein being sequenced.
    """

    def __init__(self, **kwargs):
        """
        Initializes the protein information.

        Parameters:

            product:
                The product annotation for the gene.

            dname:
                The team which identified the protein.

            gene_locus:
                The gene locus.

            sequence_prefix:
                Any prefix required when formatting the locus.
        """

        self.__kwargs = kwargs

        self.__product: str = kwargs["product"]
        self.__dname: str = kwargs["dname"]
        self.__gene_locus: str = kwargs["gene_locus"]
        self.__sequence_prefix: str = kwargs["sequence_prefix"]

        self.get_anno_info()

    def __str__(self):
        """
        Returns the string representation for the protein information.
        """

        # Check to check that all the essential information is not None
        essential_info = [self.__product, self.__dname,
                          self.__gene_locus, self.__sequence_prefix]

        if any(essential == None for essential in essential_info):
            raise errors.InsufficientProteinInfo(protein=self)

        compilation_list = []

        product_str = "product {product}".format(product=self.__product)
        compilation_list.append(product_str)

        protein_id_str = "protein_id gnl|{dname}|{gene_locus}".format(
            dname=self.__dname,
            gene_locus=self.__gene_locus
        )
        compilation_list.append(protein_id_str)

        transcript_id_str = "transcript_id gnl|{dname}|{sequence_prefix}.{gene_locus}".format(
            dname=self.__dname,
            sequence_prefix=self.__sequence_prefix,
            gene_locus=self.__gene_locus
        )

        compilation_list.append(transcript_id_str)

        # Try and get the external reference
        db_xref_str = self.get_external_reference()

        if db_xref_str is not None:
            compilation_list.append(db_xref_str)

        return '\n'.join(compilation_list)

    def __repr__(self):

        repr_ = "product:{product}\ndname:{dname}\ngene_locus:{gene_locus}\nsequence_prefix:{sequence_prefix}".format(
            product=self.__product,
            dname=self.__dname,
            gene_locus=self.__gene_locus,
            sequence_prefix=self.__sequence_prefix,
        )

        return repr_

    def get_anno_info(self):
        """
        Extracts optional information from the annotation file.
        """

        # Create a quick reference to the annotations dataframe
        anno_df = self.__kwargs["annotation_df"]

        # Retrieve the gene id
        gene_id = self.__kwargs["gene_id"]

        try:
            # Try to retrieve the optional information corresponding to this
            # protein
            anno_row = anno_df.loc[gene_id]

        except KeyError:
            # If no external reference exists a key error will be thrown. Just
            # return None if this is the case
            return

        self.extract_value(anno_row, "gene_name", add_to_dict=True)
        self.extract_value(anno_row, "gene_synonym", add_to_dict=True)

        accession = self.extract_value(
            anno_row, "accession", add_to_dict=False)

        if accession is not None:
            self.process_accession(accession)

    def extract_value(self, anno_row, value, add_to_dict=False) -> str:
        """
        Extracts a certain value from a certain row of the annotation file.

        Parameters:
            anno_row
                The row of the annotation file to extract the value.

            value:
                The value to be extracted from the row.

            add_to_dict:
                If True, the extracted value is added to self.__kwargs without
                any further processing.

        Returns:
            Returns the extracted value as a string. None otherwise.
        """

        try:
            extracted_value = anno_row[value]
        except KeyError:
            return None

        # Don't add the gene name if it has na
        if extracted_value.lower() == 'na':
            return None

        if add_to_dict:
            self.__kwargs[value] = extracted_value

        return extracted_value

    def process_accession(self, accession):
        """
        Processes the accession value from the annotation file and adds it to
        self.__kwargs.
        """

        db_xref_prefix, db_xref_id, *_ = accession.split(sep='|')

        if db_xref_prefix == "sp":
            self.__kwargs["accession"] = 'db_xref="UniProtKB/Swiss-Prot:{0}"'.format(
                db_xref_id)
        elif db_xref_prefix == "tr":
            self.__kwargs["accession"] = 'db_xref="UniProtKB/TrEMBL:{0}"'.format(
                db_xref_id)

        raise NotImplementedError(
            "Not equipped to handle db xref prefix: " + db_xref_prefix)

    def get_external_reference(self) -> str:
        """
        Attempts to find an external database reference from the annotations 
        file.

        Return:
            The appropriately formatted reference string IF a reference is 
            found, None otherwise.

        Example:
            db_xref="UniProtKB/Swiss-Prot:P12345"
        """

        # Create a quick reference to the annotations dataframe
        anno_df = self.__kwargs["annotation_df"]

        # Retrieve the gene id
        gene_id = self.__kwargs["gene_id"]

        anno_db_xref = None

        try:
            # Try to retrieve the corresponding external reference from the
            # annotation file by indexing the dataframe created using the
            # annotations file
            anno_db_xref = anno_df.loc[gene_id][1]

        except KeyError:
            # If no external reference exists a key error will be thrown. Just
            # return None if this is the case
            return None

        db_xref_prefix, db_xref_id, *_ = anno_db_xref.split(sep='|')

        if db_xref_prefix == "sp":
            return 'db_xref="UniProtKB/Swiss-Prot:{0}"'.format(db_xref_id)
        elif db_xref_prefix == "tr":
            return 'db_xref="UniProtKB/TrEMBL:{0}"'.format(db_xref_id)

        raise NotImplementedError(
            "Not equipped to handle db xref prefix: " + db_xref_prefix)


class AbstractGene(abc.ABC):
    """
    Defines an abstract gene class within a sequence.
    """

    def __init__(self, type_: str, sequence_prefix: str, gff: Gff3,
                 gene_num: int, line_start: int, line_end: int, **kwargs):
        """
        Creates a new instance of an abstract gene.

        Parameters:
            type_:
                The gene type as a string.

            gff:
                An instance of the Gff3 parser for the gff file of interest.

            gene_ID:
                The unique string identifier for the gene.

            line_start:
                The starting line of the gene sequence.

            line_end:
                The ending line of the gene sequence.
        """

        self.__type: str = type_
        self.__sequence_prefix: str = sequence_prefix
        self._gff: Gff3 = gff
        self.__gene_num: int = gene_num
        self.__line_start: int = line_start
        self.__line_end: int = line_end
        self.__kwargs = kwargs

    @property
    def gene_id(self) -> str:
        """Retrieves the gene id."""
        raise NotImplementedError

    def get_type(self) -> str:
        """Returns the segment type as a string"""

        return self.__type

    def get_lbp_list(self) -> List[Tuple[int, int]]:
        """
        Returns a list of all the location based pairs for the gene.

        Example:
            [ (5616, 26521) ]
        """

        lpb_list: list = []

        for line_num in range(self.__line_start, self.__line_end + 1):

            if self._gff.lines[line_num]['type'] == "exon":
                continue

            elif self._gff.lines[line_num]['type'] == self.__type:
                lpb_list.append(
                    (self._gff.lines[line_num]['start'], self._gff.lines[line_num]['end']))

        return lpb_list

    def get_protein_info(self) -> str:
        """
        Builds the string representation for the protein infomation.
        """

        # Make a shallow copy of the dictionary
        protein_info_dict = {**self.__kwargs}

        protein_info_dict["product"] = "PROTEIN_ANNOTATION"
        protein_info_dict["gene_num"] = self.__gene_num
        protein_info_dict["gene_locus"] = "{locus_prefix}_{gene_num:04d}".format(
            locus_prefix=protein_info_dict["locus_prefix"],
            gene_num=self.__gene_num
        )
        protein_info_dict["parent"] = self._gff.lines[self.__line_start]
        protein_info_dict["gene_id"] = self.gene_id

        simple_protein_info = ProteinInformation(
            **protein_info_dict,
            sequence_prefix=self.__sequence_prefix
        )

        return str(simple_protein_info)

    def __str__(self) -> str:
        """
        Creates a string representation for this class.
        """

        lpb_list = self.get_lbp_list()

        first_line = "{start}\t{end}\t{type_}".format(
            start=lpb_list[0][0], end=lpb_list[0][1], type_=self.__type
        )

        locs_list = []

        for start, end in lpb_list[1:]:
            locs_list.append("{start}\t{end}".format(
                start=start, end=end))

        protein_info = '\t\t\t'.join(self.get_protein_info().splitlines(True))

        if len(locs_list):
            return first_line + '\n' + '\n'.join(locs_list) + '\n' + '\t\t\t' + protein_info

        return first_line + '\n' + '\t\t\t' + protein_info


class mRNAGene(AbstractGene):
    """
    Creates a container class for a mRNA gene.
    """

    def __init__(self, gff: Gff3, gene_num: int, line_start: int, line_end: int, **kwargs):
        """
        Creates a new instace of an abstract gene.

        Arguments the same as AbstractGene documentation
        """

        super().__init__("mRNA", "mrna", gff, gene_num, line_start, line_end, **kwargs)

    @property
    def gene_id(self):
        return self._gff.lines[self._AbstractGene__line_start]["attributes"]["Name"]


class CDSGene(AbstractGene):
    """
    Creates a container class for a CDS gene.
    """

    def __init__(self, gff: Gff3, gene_num: int, line_start: int, line_end: int, **kwargs):
        """
        Creates a new instace of an abstract gene.

        Arguments the same as AbstractGene documentation
        """

        super().__init__("CDS", "mrna", gff, gene_num, line_start, line_end, **kwargs)

    @property
    def gene_id(self):
        return self._gff.lines[self._AbstractGene__line_start]["attributes"]["Parent"][0]


class GeneSequence:
    """
    Holds information about a single feature.
    """

    def __init__(self, gff: Gff3, gene_num: int,
                 line_start: int, line_end: int, **kwargs):
        """
        Creates an instance of GeneSequence.

        Parameters:
            gff:
                An instance of the Gff3 parser for the gff file of interest.

            gene_num:
                The index of the gene within the gff file.

                Example, if Bmin.gene1 is the first gene to appear in the gff
                file then its gene_num is 1.

            gene_ID:
                The unique string identifier for the gene.

            line_start:
                The starting line of the gene sequence.

            line_end:
                The ending line of the gene sequence.
        """

        self.__gff: Gff3 = gff
        self.__gene_num: int = gene_num
        self.__line_start: int = line_start
        self.__line_end: int = line_end
        self.__kwargs = kwargs

    def analyse_sequence(self) -> Dict[Tuple[int, int], str]:
        """
        Analyses the gene sequence and breaks it down into its constituent
        parts (CDS, mRNA etc).
        """

        segment_dict: Dict[Tuple[int, int], str] = {}

        current_line_num: int = self.__line_start

        # Keep moving to the next line until we find a type that is neither
        # exon nor gene

        while self.__gff.lines[current_line_num]["type"] in ("gene", "exon"):
            current_line_num += 1

        current_type = self.__gff.lines[current_line_num]["type"]

        next_type = None

        # Start decomposing the sequence
        for next_line_num in range(current_line_num + 1, self.__line_end):

            # Ignore lines with type exons
            if self.__gff.lines[next_line_num]["type"] == "exon":
                continue

            next_type = self.__gff.lines[next_line_num]["type"]

            # Check to see if we have reached a new segment
            if self.segment_condition_met(current_line_num, next_line_num):
                segment_dict[(current_line_num, next_line_num - 1)
                             ] = current_type

                current_type = next_type
                current_line_num = next_line_num

        # Don't forget to save the last segment!
        segment_dict[(current_line_num, self.__line_end)] = next_type

        return segment_dict

    def segment_condition_met(self, current_line_num: int, next_line_num: int) -> bool:
        """
        A function that checks if the next line starts a new segment within
        the gene sequence.

        Parameters:
            current_line_num:
                The starting line of the current segment that is being produced

            next_line_num:
                The next line number of the gene sequence to process.

        Return:
            Returns True if the next line defines a new segment.
        """

        # Ignore lines with type exons
        if self.__gff.lines[next_line_num]["type"] == "exon":
            return False

        # Start a new segment if the types of two consecutive lines are not
        # the same
        if self.__gff.lines[current_line_num]["type"] != \
                self.__gff.lines[next_line_num]["type"]:
            return True

        return False

    def __str__(self) -> str:
        """Creates a string representation for the gene sequence."""

        # Get the gene dictionary for this sequence
        gene_dict: Dict[Tuple[int, int], str] = self.analyse_sequence()

        # Holds all the string representations of the gene segments in a list
        segment_str_list: List[str] = []

        for type_, (gene_con_start, gene_con_end) in zip(gene_dict.values(), gene_dict.keys()):

            gene_container: AbstractGene = None

            if type_ == "CDS":
                gene_container = CDSGene(self.__gff, self.__gene_num,
                                         gene_con_start, gene_con_end, **self.__kwargs)
            elif type_ == "mRNA":
                gene_container = mRNAGene(self.__gff, self.__gene_num,
                                          gene_con_start, gene_con_end, **self.__kwargs)

            if gene_container is not None:
                segment_str_list.append(str(gene_container))

        start_end_line: str = "{start}\t{end}\tgene".format(
            start=self.__gff.lines[self.__line_start]['start'],
            end=self.__gff.lines[self.__line_start]['end']
        )

        locus_line: str = "\t\t\tlocus_tag\t{locus_prefix}_{gene_num:04d}".format(
            locus_prefix=self.__kwargs["locus_prefix"],
            gene_num=self.__gene_num
        )

        segment_str_list.insert(0, locus_line)
        segment_str_list.insert(0, start_end_line)

        return '\n'.join(segment_str_list)


class Feature:
    """
    A container for feature information.
    """

    def __init__(self, feature_name, start_gene_count, gff, line_start, line_end, **kwargs):
        """
        Initializes a new feature.
        """

        self.__feature_name: str = feature_name
        self.__start_gene_count: int = start_gene_count
        self.__end_gene_count: int = start_gene_count
        self.__gff: Gff3 = gff
        self.__line_start: int = line_start
        self.__line_end: int = line_end
        self.__kwargs: dict = kwargs

    def get__end_gene_count(self) -> int:
        return self.__end_gene_count

    def segment_condition_met(self, current_line_num, next_line_num):

        if self.__gff.lines[next_line_num]['type'] == "gene":
            return True

        return False

    def analyse(self) -> Dict[str, Tuple[int, int]]:
        """
        Creates a dictionary which contains the lines (zero indexed) within the gff file 
        which each gene corresponds to.

        Example:
            {
                "Bmin.gene1": (0,55),
                "Bmin.gene2": (56,165)
            }
        This example tells us that gene "Bmin.gene1" extends from line 0 to 
        line 55 (zero indexed) within the gff file.

        Return:
            The aforementioned dictionary.
        """

        gene_lines_dict: Dict[Tuple[int, int], str] = {}
        current_line_num: int = self.__line_start

        # Keep moving to the next line until we find a gene

        while self.__gff.lines[current_line_num]["type"] != "gene":
            current_line_num += 1

        current_gene = self.__gff.lines[current_line_num]["attributes"]["ID"]

        # Start decomposing the sequence
        for next_line_num in range(current_line_num + 1, self.__line_end):

            # Check to see if we have reached a new segment
            if self.segment_condition_met(current_line_num, next_line_num):
                gene_lines_dict[(current_line_num, next_line_num - 1)
                                ] = current_gene

                current_line_num = next_line_num
                current_gene = self.__gff.lines[next_line_num]["attributes"]["ID"]

        # Don't forget to save the last segment!
        gene_lines_dict[(current_line_num, self.__line_end)] = current_gene

        return gene_lines_dict

    def __str__(self):
        """
        Creates a string representation for the Feature.
        """

        # Create the feature header
        feature_header: str = ">Feature {feature_name}".format(
            feature_name=self.__feature_name)

        # Make a list of components that will eventually
        comp_list = [feature_header]

        gene_lines_dict: Dict[Tuple[int, int], str] = self.analyse()

        for gene_num, (gene_start, gene_end) in enumerate(gene_lines_dict.keys(),
                                                          start=self.__start_gene_count):

            new_gene = GeneSequence(
                self.__gff, gene_num, gene_start, gene_end, **self.__kwargs)
            comp_list.append(str(new_gene))

        self.__end_gene_count = gene_num

        return '\n'.join(comp_list)


class FlatFileCreator:
    """
    A class that manages the overall creation of the final flat file.
    """

    __gff_path = path.FileExtPath("gff")
    __annotation_path = path.FilePath()

    def __init__(self, locus_prefix: str, dname: str, gff_path: str,
                 annotation_path: str, *, output_path: str = None, anno_delim: str = '\t'):
        """
        Constructs a new Flat File Constructor.

        Parameters:
            locus_prefix:
                The prefix of each locus.

            dname:
                The teamname that processed the data.

            gff_path:
                A file path to the gff file.

            output_path:
                A file path to output the contents of the flatfile.

            anno_delim:
                The delimiter for the annotation file. The default is a tab.
        """

        if anno_delim == 'None':
            anno_delim = None

        else:
            for old, new in [('\\n', '\n'), ('\\t', '\t'), ('\\r', '\r')]:
                anno_delim = anno_delim.replace(old, new)

        self.__gff_path: str = gff_path
        self.__anno_delim: str = anno_delim
        self.__annotation_path: str = annotation_path
        self.__output_path: str = output_path

        # Open the gff file and read it into a dict
        self.__gff: Gff3 = Gff3(gff_file=self.__gff_path)

        # Open the annotations path as a pandas dataframe.
        # Set sep=None so that the Python parsing engine can automatically
        # detect the separator.
        self.__annotation_df = self.open_anno_file()

        # Start creating kwarg dict for protein information
        self.__kwarg = {
            "dname": dname,
            "locus_prefix": locus_prefix,
            "annotation_df": self.__annotation_df
        }

    def create_flatfile(self):
        """
        Creates a flatfile based on the input gff file.
        """

        if self.__output_path is None:
            print(str(self), flush=True)
            return

        with open(self.__output_path, "w") as file_out:
            print(str(self), file=file_out, flush=True)

        return

    def open_anno_file(self):
        """
        Opens the annotations file.
        """

        # Try opening the annotations file using the given delimiter using the
        # given delimiter using the c-engine. We will only use the python engine
        # if specified the delimiter is None.

        if self.__anno_delim is None:
            return pd.read_csv(
                self.__annotation_path, engine='python', sep=None, header=None, index_col=0)

        anno_df = pd.read_csv(
            self.__annotation_path, engine='c', sep=self.__anno_delim, header=None, index_col=0)

        _, df_cols = anno_df.shape

        if df_cols != 0:
            return anno_df

        # If we only have one column then the delimiter we are using probably is
        # not correct. We should have at least two columns.

        warning_message = 'Could not separate annotation dataframe' + \
            'with --delimiter={0!r}. Switching to python engine parser.'
        warning_message = warning_message.format(self.__anno_delim)

        warnings.warn(warning_message, RuntimeWarning)

        return pd.read_csv(
            self.__annotation_path, engine='python', sep=None, header=None, index_col=0)

    def segment_condition_met(self, current_line_num, next_line_num):

        current_feature_name = self.__gff.lines[current_line_num]['seqid']
        next_feature_name = self.__gff.lines[next_line_num]['seqid']

        if current_feature_name != next_feature_name:
            return True

        return False

    def analyse(self) -> Dict[Tuple[int, int], str]:
        """
        Creates a dictionary which contains the start and end lines of each new
        feature.

        Example:
            {
                (0, 901): "Bmin.scaffold2",
                (902, 1259): "Bmin.scaffold2",
                (1260, 1951): "Bmin.scaffold3",
            }
        This example tells us that gene "Bmin.scaffold2" extends from line 0 to 
        line 901 (zero indexed) within the gff file.

        Return:
            The aforementioned dictionary.
        """

        feat_lines_dict: Dict[str, Tuple[int, int]] = {}

        current_feature, current_line_num = self.__gff.lines[0]['seqid'], 0

        for next_line_num in range(len(self.__gff.lines)):

            if self.segment_condition_met(current_line_num, next_line_num):
                feat_lines_dict[(current_line_num,
                                 next_line_num - 1)] = current_feature

                current_line_num = next_line_num
                current_feature = self.__gff.lines[next_line_num]['seqid']

        feat_lines_dict[(current_line_num,
                         len(self.__gff.lines) - 1)] = current_feature

        return feat_lines_dict

    def __str__(self):
        """
        Creates the string output for the flatfile.
        """

        # Gets the lines within the gff file which each gene belongs to
        feat: Dict[str, Tuple[int, int]] = self.analyse()

        # Hold all the string representations of genes in a list
        feat_str_list: list = []

        start_gene_count = 1

        for feat_name, (start, end) in zip(feat.values(), feat.keys()):

            tmp_feat = Feature(feat_name, start_gene_count, self.__gff,
                               start, end, **self.__kwarg)
            feat_str_list.append(str(tmp_feat))

            start_gene_count = tmp_feat.get__end_gene_count() + 1

        return '\n'.join(feat_str_list)
