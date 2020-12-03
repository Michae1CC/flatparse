#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import abc
import os


class AutoFilePath:

    __counter = 0

    def __init__(self, target_dir: str = os.getcwd()):
        cls = self.__class__
        prefix = cls.__name__
        index = cls.__counter
        self.storage_name = '_{}#{}'.format(prefix, index)
        cls.__counter += 1

        self.target_dir = target_dir

    def __get__(self, instance, owner):
        """
        Retrieves the file path string from a manager instance.
        """

        if instance.__dict__[self.storage_name] is None:
            raise ValueError(
                f"The path requested has not been set to a legal value and is still 'None'.")

        return instance.__dict__[self.storage_name]

    def __set__(self, instance, value):
        """
        Returns the stored key value of the __dict__ attribute without any checks.
        """
        instance.__dict__[self.storage_name] = value


class ValidatedFilePath(abc.ABC, AutoFilePath):

    def __set__(self, instance, value):
        """
        Returns the stored key value of the __dict__ attribute with an
        additional validation check.
        """
        value = self.validate(instance, value)
        super().__set__(instance, value)

    @abc.abstractmethod
    def validate(self, instance, value):
        """Returns a validated value or raises an error"""


class FilePath(ValidatedFilePath):
    """The store value must be a valid file path"""

    def validate(self, instance, value):
        """
        Sets a new file path to a file path instance.
        """

        if value is None:
            # The path has not been set yet
            return None
        elif os.path.isfile(value):
            # Get the absolute path to the directory
            return os.path.abspath(value)
        elif os.path.isfile(os.path.join(self.target_dir, value)):
            # Try adding the target directory to the path name
            return os.path.abspath(
                os.path.join(self.target_dir, value))

        # This is not a valid path, throw an exception
        raise OSError(
            f"The specified filepath file could not be found.\n{value}")


class FileExtPath(FilePath):
    """
    The store value must be a valid file path with the correct 
    extension type.
    """

    def __init__(self, file_extension: str, target_dir: str = os.getcwd()):
        """
        Parameters:
            file_extension:
                A singleton or iterable of file extensions the target_dir must use
                as a suffix.
        """
        super().__init__(target_dir=target_dir)

        # Save the given file extension
        self.file_extension = list(file_extension)

    def validate(self, instance, value):

        if value is not None:
            basename = os.path.basename(value)

            # Check if the the basename end withs any of the extensions within
            # self.file_extension
            if not any(basename.endswith(ext) for ext in self.file_extension):
                raise OSError(f"The base name of the specified file '{basename}' "
                              f"does not have the correct extension type '{self.file_extension}'.")

        # Use the validation from the FilePath
        return super().validate(instance, value)


class FolderPath(ValidatedFilePath):
    """The store value must be a valid file path"""

    def validate(self, instance, value):
        """
        Sets a new folder path to a folder path instance.
        """

        if value is None:
            # The path has not been set yet
            return None
        elif os.path.isdir(value):
            # Get the absolute path to the directory
            return os.path.abspath(value)
        elif os.path.isdir(os.path.join(self.target_dir, value)):
            # Try adding the target directory to the path name
            return os.path.abspath(
                os.path.join(self.target_dir, value))

        # This is not a valid path, throw an exception
        raise OSError(
            f"The specified folder path could not be found.\n{value}")
