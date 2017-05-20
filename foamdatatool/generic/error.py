#!/usr/bin/python3
"""
Name:           error.py
Last Modified:  05/20/2017
Python Version: 3.4.3

Description:
    Module providing exception objects specific to parsing OpenFOAM data.

Contents:
    class Error:   
        A generic class for error handling
    class PathError:
        Class for handling file path errors
    class DataTypeError:
        Class for handling data type errors
    class PatternFailure:   
        Class for handling regex failures

Usage:
    These classes are used to extend the functionality of the standard Exception
    object.  All exception handling is performed in the usual manner.  However,
    using a description string is encouraged.

Notes:
"""

__all__ = ["Error", "PathError", "DataTypeError", "PatternFailure"]

# Custom exception base class
class Error(Exception):
    pass

# Custom exception sub-class for invalid file paths
class PathError(Error):
    def __str__(self):
        return str("The following path is invalid:"
                   " {}".format(Error.__str__(self)))

# Custom exception sub-class for data determination failure
class DataTypeError(Error):
    def __str__(self):
        return str("Data type determination failed in function:"
                   " {}".format(Error.__str__(self)))

# Custom exception sub-class for regular expression failures
class PatternFailure(Error):
    def __string(self):
        return str("The following pattern failed to produce expected"
                   " results: {}".format(Error.__str__(self)))
