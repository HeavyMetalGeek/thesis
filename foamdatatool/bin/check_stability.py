#!/usr/bin/python3
"""
Name:           check_stability.py
Last Modified:  05/20/2017
Python Version: 3.4.3

Description:
    Module used for calculating and displaying the autocorrelation of velocity
    components at a probed point as the time lag increases.  Instantaneous
    values are graphed for comparison.

Contents:
    class PointData :   A generic class for encapsulating data from an OpenFOAM
                        point probe data file
    func main       :   The main function.

Usage:
    main(f, start=0, end=-1, freq=10, lagmax=2)
        f       :   File containing the probe data
        start   :   Index of first data point
        end     :   Index of last data point
        freq    :   Sampling rate in hertz
        lagmax  :   Maximum time lag

Notes:
    If called as __main__, computations will be conducted using the parameters
    set within the if-statement at the bottom of the module
"""

import sys
import os
sys.path.insert(0, os.path.abspath('../..'))
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

class PointData:
    def __init__(self):
        pass

def main(f, start=0, end=-1, freq=10, lagmax=2):
    # Extract all time and point data
    data = []
    with open(f, 'r') as fname:
        for i in fname.readlines():
            temp = i.split()
            if(temp[0] == '#'):
                continue
            for j in range(len(temp)):
                temp[j] = temp[j].strip(' ')
                temp[j] = temp[j].strip('(')
                temp[j] = temp[j].strip(')')
                temp[j] = float(temp[j])
            data.append(temp)
    data = np.array(data)
    
    if(start != 0):
        start = np.where(data[:,0] == start)[0][0]
    if(end != -1):
        end = np.where(data[:,0] == end-lagmax*freq)[0][0]
    
    ptdata = PointData()
    ptdata.data = data
    ptdata.freq = freq
    
    ptdata.t = ptdata.data[start:end,0]
    ptdata.Ux = ptdata.data[start:end,1]
    ptdata.Uy = ptdata.data[start:end,2]
    ptdata.Uz = ptdata.data[start:end,3]
    
    Ruu = []
    Rvv = []
    Rww = []

    tlag = []
    
    ptdata.t0 = ptdata.t[0]
    ptdata.t1 = ptdata.t[-lagmax*freq]
    
    for i in range(0, lagmax*freq + 1):
        tlag.append(i / freq)
        Rxx = 0
        Ryy = 0
        Rzz = 0
        for j in range(0,len(ptdata.t)-lagmax*freq):
            Rxx += 1/(ptdata.t1-ptdata.t0) * ptdata.Ux[j] * ptdata.Ux[j+i]
            Ryy += 1/(ptdata.t1-ptdata.t0) * ptdata.Uy[j] * ptdata.Uy[j+i]
            Rzz += 1/(ptdata.t1-ptdata.t0) * ptdata.Uz[j] * ptdata.Uz[j+i]
        Ruu.append(Rxx)
        Rvv.append(Ryy)
        Rww.append(Rzz)
    
    Ruu = np.array(Ruu)
    Rvv = np.array(Rvv)
    Rww = np.array(Rww)
    Cuu = Ruu / Ruu[0]
    Cvv = Rvv / Rvv[0]
    Cww = Rww / Rww[0]
        
    #props = dict(boxstyle='round', facecolor='wheat', alpha=1)
    fig1 = plt.figure(1, figsize=(20,10))
    #ax1 = fig1.add_axes([0.1,0.1,0.8,0.8])
    gs = GridSpec(2,1,top=0.95, bottom=0.05, left=0.05, right=0.95, hspace=0.25, height_ratios=[3,1])
    ax1 = fig1.add_subplot(gs[0])
    ax1
    ax1.set_title(r'Auto-correlation from t = {} to {}'.format(
            ptdata.t0, ptdata.t1))
    ax1.set_xlabel(r'Time Lag [s]')
    ax1.set_ylabel(r'Auto-correlation')
    ax1.set_ylim(ymax=1.005, ymin=0.975)
    '''
    ax1.text(0.05, 0.95, 
             r'SMG $U_x = {:5.3f}$'.format(np.mean(oData[:,1])) + 
             r'1EQ $U_x = {:5.3f}$'.format(np.mean(nData[:,1])),
            horizontalalignment='left',
            verticalalignment='top',
            size=12,
            transform=ax1.transAxes,
            bbox=props)
    '''
    ax1.plot(tlag, Cuu, marker='.', label=r'$u$')
    ax1.plot(tlag, Cvv, marker='+', label=r'$v$')
    ax1.plot(tlag, Cww, marker='x', label=r'$w$')
    ax1.legend()
    
    ax2 = fig1.add_subplot(gs[1])
    ax2.set_title(r'Instantaneous Velocity at Mesh Center')
    ax2.set_xlabel(r'Time [s]')
    ax2.set_ylabel(r'Velocity [m/s]')
    ax2.plot(ptdata.t, ptdata.Ux, label=r'$u^\prime$')
    ax2.plot(ptdata.t, ptdata.Uy, ls='--', label=r'$v^\prime$')
    ax2.plot(ptdata.t, ptdata.Uz, ls=':', label=r'$w^\prime$')
    #ax2.axhline(np.mean(ptdata.Ux), ls='-.', c='black')
    ax2.legend()
    
    plt.show()

if(__name__ == '__main__'):
    ROOT = '/media/fatirishman53/96A3-CF1A/data'
    myfile = ROOT+'/sm_nW/nutU/postProcessing/centerProbe/0/U'
    main(myfile, start=4000, freq=10, lagmax=2)