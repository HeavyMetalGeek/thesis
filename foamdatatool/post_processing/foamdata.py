#!/usr/bin/python3
"""
Name:           FoamData.py
Last Modified:  05/20/2017
Python Version: 3.4.3

Description:
    Module containing classes for parsing, storing, and processing OpenFOAM
    data.

Contents:
    class FoamData:
        A class for parsing, storing, and processing OpenFOAM data
    class VelocityData:
        Sub-class of Foam Data specifically for velocity

Usage:
    FoamData(time, valname, casepath=os.getcwd())    
        instance.data:   
            Array of retrieved data
        instance.get_data():
            Retrieve data from valname file
        instance():
            Calls instance.get_data()
        instance.get_unique():   
            Creates array of unique values
        instance.unique:
            List of unique values
        print(instance):
            Print path and length of data
    
    VelocityData(time, casepath=os.getcwd())
        instance.get_data(component=None):
            Retreives velocity data.  Specifying a component will narrow data
            to that specific component (x=1, y=2, z=3)
        instance():
            Calls instance.get_data()

Notes:
"""

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))
import re
import numpy as np
from foamdatatool.generic.error import Error
from foamdatatool.generic.error import PathError
from foamdatatool.generic.error import DataTypeError
from foamdatatool.generic.error import PatternFailure


__all__ = ["FoamData", "VelocityData"]

# Base class for OpenFOAM data retrieval
class FoamData:
    def __init__(self, time, valname, casepath=os.getcwd()):
        self.re_foam_time = re.compile(r'^[0-9]+\.*[0-9]*$')
        self.re_nulist_scalar = re.compile(
            r'(?:internalField[ \n]+nonuniform List\<scalar\>'
            r'[ \n]*[0-9]+[ \n]*\()(.*?\)[ \n]*;)', re.M|re.S)
        self.re_nulist_vector = re.compile(
            r'(?:internalField\s+nonuniform List\<vector\>'
            r'[ \n]*[0-9]+[ \n]*\()(.*?\)[ \n]*;)', re.M|re.S)
        self.re_scalar = re.compile(r'(^[-]?[0-9]+.*[0-9]*$)', re.M)
        self.re_vector = re.compile(
            r'^(?:\(\s*)([+-]?[0-9]+\.*[0-9e-]*\s[+-]?[0-9]+\.*'
            r'[0-9e-]*\s[+-]?[0-9]+\.*[0-9e-]*)(?:\s*\))$', re.M)
        self.casepath = str(casepath)
        self.time = str(time)
        self.valname = str(valname)
        self.datatype = None
        self.datapath = (self.casepath 
                         + '/' + self.time 
                         + '/' + self.valname)
        self.length = 0
    
    # BULLET: Consider added functionality
    def __call__(self):
        self.get_data()
    
    # Produces a printout of the data path and actual data   
    def __str__(self):
        return str(str(self.data) 
                    + '\n\n  FROM: ' + self.datapath
                    + '\nLENGTH: ' + str(self.length))
        
    def get_data(self):
        try:
            if(not os.path.exists(self.datapath)): 
                raise PathError(self.datapath)
            with open(self.datapath, 'r') as infile:
                full_string = infile.read()
                is_vector = self.re_nulist_vector.search(full_string)
                is_scalar = self.re_nulist_scalar.search(full_string)
                if(is_vector == is_scalar):
                    raise DataTypeError('FoamData.get_data()')
                elif(is_vector):
                    self.datatype = 'non-uniform list vector'
                    print("Extracting vectors in {}"
                          .format(self.datapath))
                    all_vals = self.re_nulist_vector.findall(
                               full_string)
                    if(len(all_vals) == 0): 
                        raise PatternFailure('all_vals' +
                            self.re_nulist_vector.pattern)
                    sep_vals = self.re_vector.findall(all_vals[0])
                    if(len(sep_vals) == 0):
                        raise PatternFailure('sep_vals' +
                            self.re_vector.pattern)
                    for i in range(len(sep_vals)):
                        sep_vals[i] = re.sub(' {1}', ',', sep_vals[i])
                        sep_vals[i] = [float(k) for k in 
                                       sep_vals[i].split(',')]
                    self.data = np.array(sep_vals)
                elif(is_scalar):
                    self.datatype = 'non-uniform list scalar'
                    print("Extracting scalars in {}"
                          .format(self.datapath))
                    all_vals = self.re_nulist_scalar.findall(
                               full_string)
                    if(len(all_vals) == 0): 
                        raise PatternFailure(
                            self.re_nulist_scalar.pattern)
                    sep_vals = self.re_scalar.findall(all_vals[0])
                    sep_vals = [float(k) for k in sep_vals]
                    if(len(sep_vals) == 0):
                        raise PatternFailure(
                            self.re_scalar.pattern)
                    self.data = np.array(sep_vals)
                else:
                    raise DataTypeError('FoamData.get_data()')
                self.length = len(self.data)    
        except DataTypeError as err:
            print(err)
        except PatternFailure as err:
            print(err)
        except PathError as err:
            print(err)
                

    def get_unique(self):
        self.unique = []
        for i in range(len(self.data)):
            if not(self.data[i] in self.unique):
                self.unique.append(self.data[i])
        self.unique = np.array(self.unique)

# Sub-class for specifically handling velocity values 
class VelocityData(FoamData):
    def __init__(self, time, casepath=os.getcwd(), 
                 valname='U', component=None):
        FoamData.__init__(self, time, valname, casepath)
        self.component = component
        
    def __call__(self):
        self.get_data()
    
    # Method for retrieving velocity data from OpenFOAM files
    def get_data(self, component=None):
        try:
            if(not os.path.exists(self.datapath)): 
                raise PathError(self.datapath)
            with open(self.datapath, 'r') as infile:
                full_string = infile.read()
                is_vector = self.re_nulist_vector.search(full_string)
                if(is_vector):
                    self.datatype = 'non-uniform list vector'
                    print("Extracting vectors in {}"
                          .format(self.datapath))
                    all_vals = self.re_nulist_vector.findall(
                               full_string)
                    if(len(all_vals) == 0): 
                        raise PatternFailure('all_vals' +
                            self.re_nulist_vector.pattern)
                    sep_vals = self.re_vector.findall(all_vals[0])
                    if(len(sep_vals) == 0):
                        raise PatternFailure('sep_vals' +
                            self.re_vector.pattern)
                    for i in range(len(sep_vals)):
                        sep_vals[i] = re.sub(' {1}', ',', sep_vals[i])
                        sep_vals[i] = [float(k) for k in 
                                       sep_vals[i].split(',')]
                        # Keep only specified component
                        if(component in [1,2,3]):
                            sep_vals[i] = sep_vals[i][component-1]
                    self.data = np.array(sep_vals)
                    self.length = len(self.data)
                else:
                    raise Error('NO VELOCITY DATA FOUND')
        except Error as err:
            print(err)
