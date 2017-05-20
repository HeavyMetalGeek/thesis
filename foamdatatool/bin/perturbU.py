#!/usr/bin/python3
"""
Name:           perturbU.py
Last Modified:  05/20/2017
Python Version: 3.4.3

Description:
    Module used to test the functionality of the perturbU utility in OpenFOAM.

Contents:
    class FoamVector:
        Class representing an OpenFOAM vector quantity
    class FoamScalar:
        Class representing an OpenFOAM scalar quantity
    func mag:
        Returns the magnitude of a vector
    func getCellCentres:
        Returns the cell center coordinates


Usage:
    FoamVector(x, y, z)
        instance.x:
            X-component of the vector
        instance.y:
            Y-component of the vector
        instance.z:
            Z-component of the vector
        instance.value:
            Returns numpy array of vector components

    FoamScalar(val)
        instance.val:
            Scalar value
        instance.value:   
            Returns numpy array of scalar value
        
    mag(f, start=0, end=-1, freq=10, lagmax=2)
        vector:   
            Numpy array representing a vector quantity.  Must be a 1x3 
            array to make sense in this context

    getCellCentres(f)
        f:
            Case file containing cell center data
            *Note that this requires that writeCellCentres has been run already, 
            producing the ccx, ccy, and ccz files

Notes:
    MODULE CLEANING NEEDED!
    *Encapsulate processes into functions
    *Modify for more generalized use
"""

import sys
import os
sys.path.insert(0, os.path.abspath('../..'))
import logging
import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt
from foamdatatool.post_processing.foamdata import FoamData
from foamdatatool.generic.file_strings import of_header
from foamdatatool.generic.error import Error
from matplotlib.gridspec import GridSpec

def mag(vector):
    return np.linalg.norm(vector)

def getCellCentres(f):
    x = FoamData(time=0, valname='ccx', casepath=f)
    y = FoamData(time=0, valname='ccy', casepath=f)
    z = FoamData(time=0, valname='ccz', casepath=f)
    x.get_data()
    y.get_data()
    z.get_data()
    return np.column_stack((x.data,y.data,z.data))

class FoamVector:
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        
    def value(self):
        return np.array([self.x, self.y, self.z])

class FoamScalar:
    def __init__(self, val):
        self.val = val
        
    def value(self):
        return np.array([self.val])

#ROOT = '/home/fatirishman53/OpenFOAM/fatirishman53-3.0.x/run/nist_p'
ROOT = '/home/fatirishman53/OpenFOAM/fatirishman53-4.1/run/p_nist'

yWD = FoamData(time=0, valname='ccy', casepath=ROOT)
yWD.get_data()
yWD.get_unique()
yw = yWD.data

h = np.max(yw)

yDir = np.full((len(yw.data),3), [0., -1., 0.], int)

U = np.full((len(yw.data),3), 0, float)

nu = FoamScalar(1.455e-5)
Ubar = FoamVector(27,0,0)
Retau = FoamScalar(811683.85)

xDir = Ubar.value() / mag(Ubar.value())

print("Re(tau) = ", Retau.value())
utau = Retau.value() * nu.value() / h 
print("u(tau) = ", utau)

duplus = Ubar.value()[0] * 0.25 / utau
betaPlus = 2 * np.pi * (1/200)
sigma = 0.00055
alphaPlus = 2 * np.pi * (1/500)
epsilon = Ubar.value()[0]/200

print("duplus = ", duplus)
print("betaPlus = ", betaPlus)
print("sigma = ", sigma)
print("alphaPlus = ", alphaPlus)
print("epsilon = ", epsilon)

ld = input("Load dumped centres data? ")
f = ROOT+'/centres.dump'
if(ld in ['y', 'Y', 'yes', 'Yes', 'YES']):    
    if not(os.path.exists(f)):
        raise Error("File not found: {}".format(f))
        sys.exit(0)
    with open(f, 'rb') as fname:
        centres = pickle.load(fname)
else:
    centres = getCellCentres(ROOT)
    with open(f, 'wb') as fname:
        pickle.dump(centres, fname)

perturbation = np.random.RandomState(1234567)

ld = input("Load dumped U_p data? ")
f = ROOT+'/U_p.dump'
strm_streaks = []
span_streaks = []
if(ld in ['y', 'Y', 'yes', 'Yes', 'YES']):
    if not(os.path.exists(f)):
        raise Error("File not found: {}".format(f))
        sys.exit(0)
    with open(f, 'rb') as fname:
        U = pickle.load(fname)
else:
    print("Computing U values")
    for i in range(len(centres)):
        deviation = 1 + 0.1*perturbation.normal()
        cCentre = centres[i]
        
        zDir = np.cross(xDir, yDir[i])
        zDir = zDir / np.linalg.norm(zDir)
        
        zplus = np.dot(cCentre, zDir) * Retau.value() / h
        yplus = yw[i] * Retau.value() / h;
        xplus = np.dot(cCentre, xDir) * Retau.value() / h
        
        # Laminar parabolic profile
        U[i] = 3*Ubar.value() * (yw[i]/h - 0.5*np.square(yw[i]/h))
        
        # Streak streamwise velocity
        U[i] += xDir * (utau*duplus/2) \
            * (yplus/40) \
            * np.exp(-sigma*np.square(yplus) + 0.5) \
            * np.cos(betaPlus*zplus) \
            * deviation
        strm_streaks.append(xDir * (utau*duplus/2)
            * (yplus/40)
            * np.exp(-sigma*np.square(yplus) + 0.5)
            * np.cos(betaPlus*zplus)
            * deviation)
        
        ex = True
        if(utau == 0):
            print("utau = ", utau)
            ex = True
        if(np.all(np.equal(xDir, np.array([0,0,0])))):
            print("xDir = ", xDir)
            ex = True
        if(np.all(np.equal(yplus, np.array([0,0,0])))):
            print("yplus = ", yplus)
            ex = True
        if(np.all(np.equal(zplus, np.array([0,0,0])))):
            print("zplus = ", zplus)
            ex = True
        if(deviation == 0):
            print("deviation = ", deviation)
            ex = True
        if(ex):
            print("utau = ", utau)
            print("xDir = ", xDir)
            print("yplus = ", yplus)
            print("zplus = ", zplus)
            print("term1 = ", xDir * (utau*duplus/2))
            print("term2 = ", yplus/40)
            print("term3 = ", np.exp(-sigma*np.square(yplus) + 0.5))
            print("term4 = ", np.cos(betaPlus*zplus))
            print("term5 = ", deviation)
            sys.exit(0)
        # Streak spanwise perturbation
        U[i] += zDir * epsilon \
            * np.sin(alphaPlus*xplus) \
            * yplus \
            * np.exp(-sigma*np.square(yplus)) \
            * deviation
        span_streaks.append(zDir * epsilon
            * np.sin(alphaPlus*xplus)
            * yplus
            * np.exp(-sigma*np.square(yplus))
            * deviation)
        pct = i/len(centres) * 100
        if(pct % 10 == 0):
            print("Computing... {:5.1f}%\r".format(pct), end='')
    print("Calculation Complete")
    print("Writing Data")
    strm_streaks = np.array(strm_streaks)
    span_streaks = np.array(span_streaks)
    with open(ROOT+'/streak_data', 'w') as fname:
        for i in range(1000):
            fname.write("{}\t\t{}\n".format(strm_streaks[i], span_streaks[i]))
    print("Creating dump file", f)
    with open(f, 'wb') as fname:
        pickle.dump(U, fname)

print("Writing to file", ROOT+'/0/U_p')
with open(ROOT+'/0/U_p', 'w') as fname:
    fname.write(of_header)
    fname.write('\n\n')
    fname.write("dimensions\t[0 1 -1 0 0 0 0];")
    fname.write('\n\n')
    fname.write("internalField\tnonuniform List<vector>")
    fname.write('\n')
    fname.write(str(len(U)))
    fname.write('\n')
    fname.write('(\n')
    for i in range(len(U)):
        fname.write("({} {} {})".format(U[i][0], U[i][1], U[i][2]))
        fname.write('\n')
    fname.write(')\n;\n\n')
    fname.write('boundaryField\n')
    fname.write('{\n')
    fname.write('\tbottomWall\n')
    fname.write('\t{\n\t\ttype\t\tfixedValue;\n')
    fname.write('\t\tvalue\t\tuniform (0 0 0);\n\t}\n\n')
    fname.write('\ttopWall\n')
    fname.write('\t{\n\t\ttype\t\tslip;\n\t}\n\n')
    fname.write('\tsides1_half0\n')
    fname.write('\t{\n\t\ttype\t\tcyclic;\n\t}\n\n')
    fname.write('\tsides1_half1\n')
    fname.write('\t{\n\t\ttype\t\tcyclic;\n\t}\n\n')
    fname.write('\tinout1_half0\n')
    fname.write('\t{\n\t\ttype\t\tcyclic;\n\t}\n\n')
    fname.write('\tinout1_half1\n')
    fname.write('\t{\n\t\ttype\t\tcyclic;\n\t}\n}')
    
'''
print("Ploting vectors")
from mpl_toolkits.mplot3d import axes3d

fig = plt.figure()
ax = fig.gca(projection='3d')

x,y,z = np.meshgrid(centres[:,0],
                    centres[:,1],
                    centres[:,2])

ax.quiver(x,y,z,U[:,0],U[:,1],U[:,2], length=0.1)
plt.show()
'''    