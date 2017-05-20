#!/usr/bin/python3
"""
Name:           FoamData.py
Last Modified:  05/20/2017
Python Version: 3.4.3

Description:
    Module containing a classes for parsing, processing, and storing OpenFOAM
    probe data.

Contents:
    class ProbeData:
        Class for parsing, processing, and storing OpenFOAM probed velocity data

Usage:
    ProbeData(ufile=None, ptfile=None)   
        instance.ufile:
            Path of file for velocity values
        instance.ptfile
            Path of file for point coordinates
        instance.dumpData(f,data):
            Pickle "data" and save to "f"
        instance.loadData(f):
            Load pickled data from "f" and return it
        instance.getVectors(f,probes=False):
            Processes data formatted as OpenFOAM vector data (probe vector 
            format if probes=True, non-uniform vector format if probes=False)     
        instance.processPoints(pts_file=None):
            Processes point coordinates from the path given by pts_file or the
            ptfile property.  
            *Note that if a path is provided by pts_file, it will be processed
            and the class property `ptfile` will be set to its value.
        instance.processVelocity(U_file=None):
            Processes velocity data from the path given by the U_file or ufile
            property.  
            *Note* that if a path is provided by U_file, it will be processed
            and the class property ufile will be set to its value.

Notes:
"""

import re
import logging
import numpy as np
import pickle

__all__ = ["ProbeData"]

class ProbeData:
    def __init__(self, ufile=None, ptfile=None):
        self.ufile = ufile
        self.ptfile = ptfile
        self.pts = None
        self.U = None
        if(self.ufile != None):
            self.processVelocity(ufile)
        if(self.ptfile != None):
            self.processPoints(ptfile)
        
    def __call__(self):
        pass

    def __str__(self):
        pass

    def dumpData(self, f, data):
        with open(f, 'wb') as fname:
            pickle.dump(data, fname)
    
    def loadData(self, f):
        with open(f, 'rb') as fname:
            src = pickle.load(fname)
        return src
    
    def getVectors(self, f, probes=False):
        vectors = ''
        if(probes):
            logging.debug("FORMAT: probe")
            vectors = re.compile(r'^\s+[0-9]+\.*[0-9]*\s+\(.+\)',re.S)
        else:
            logging.debug("FORMAT: non-uniform list")
            vectors = re.compile(r'^\(.+\)')
        vals = []
        with open(f, 'r') as fname:
            for i in fname.readlines():
                if not(vectors.search(i)):
                    #logging.debug("Skipping: {}".format(i.strip('\n')))
                    continue
                tempArray = i.split()
                for j in range(len(tempArray)):
                    tempArray[j] = tempArray[j].strip(' ')
                    tempArray[j] = tempArray[j].strip('(')
                    tempArray[j] = tempArray[j].strip(')')
                    tempArray[j] = float(tempArray[j])
                vals.append(tempArray)
            logging.debug("RETURNING: {} values".format(len(vals)))
            return np.array(vals)

    def processPoints(self, pts_file=None):
        if(self.ptfile == None and pts_file != None):
            self.ptfile = pts_file
        if(self.ptfile != None):
            logging.debug("PARSING POINTS: {}".format(self.ptfile))
            # Point array is point x component
            self.pts = self.getVectors(self.ptfile)
        else:
            print("No point file path to process...")
            return

    def processVelocity(self, U_file=None):
        if(self.ufile == None and U_file != None):
            self.ufile = U_file
        if(self.ufile != None):
            logging.debug("PARSING U: {}".format(self.ufile))
            self.U = self.getVectors(self.ufile, probes=True)
            # Convert velocity data into 3D array
            # Velocity array is time x point x component
            self.U = self.U[:,1:].reshape(len(self.U), -1, 3)
        else: 
            print("No velocity file path to process...")
            return