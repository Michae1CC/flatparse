#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''


class InsufficientProteinInfo(AttributeError):
    """A error raised when the protein information has not been provided with
    the minimal amount of information."""

    def __init__(self, message: str = None, protein=None):

        super_message = message

        if message is None:
            super_message = "None-type values found in ProteinInformation"

        if protein is not None:
            super_message += "\n\nEssential protein values for reference:\n" + \
                repr(protein)

        super().__init__(message)
