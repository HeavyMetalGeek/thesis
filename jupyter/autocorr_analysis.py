import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
from velocity_stats import VelocityData

def main(show=True, write_path=None):
    fine_src = './postProcessing/center_probe/290.5/U'
    cors_src = './postProcessing_60/center_probe/50/U'
    ac_plots(fine_src, show=True, write_path=None)
    mean_plots2(fine_src, window=4000, show=True, write_path=None)

def ac_plots(src, show=True, write_path=None):
    fine_src = src
    VelocityData.plot_config()
    fig, ax = plt.subplots(1)
    #x_axis_vals = np.arange(fine_800.N)*fine_800.dt
    lags = [800, 400, 200, 80, 40]
    colors = ['b', 'r', 'g', 'c', 'k']
    for i, c in zip(lags, colors):
        fine = VelocityData(fine_src)
        fine.autocorr(summation=False, N=i, lag=i)
        x_axis_vals = np.arange(fine.Rxx.size)*fine.dt
        try:
            ax.plot(x_axis_vals, fine.Rxx, color=c, label=fine.N)
        except ValueError:
            self.debug(0, "Plot axis value mismatch:")
            self.debug(0, "N={}".format(len(x_axis_vals)))
            self.debug(0, "Rxx={}".format(len(fine.Rxx)))
            return
    ax.set_title("Velocity Autocorrelation")
    ax.set_xlabel("Time Lag (s)")
    ax.set_ylabel("$R_{ii}$")
    ax.set_ylim(-1,1)
    ax.legend()
    # Write to file if write_path given
    if (write_path != None):
        folder = os.path.dirname(write_path)
        if not (os.path.exists(folder)):
            self.debug(0, "Graph could not be written to file.")
            self.debug(0, "{} does not exist".format(folder))
        else:
            fig.savefig(write_path)
    if (show):
        plt.show()

def mean_plots(src, window=80, show=True, write_path=None):
    fine_src = src
    VelocityData.plot_config()
    fig, ax = plt.subplots(1)
    
    fine = VelocityData(fine_src)
    fine.resize_data(-8000)
    fine.m_avg = pd.rolling_mean(fine.up, window=window)
    fine.m_avg = fine.m_avg[window:]
    x_axis_vals = np.arange(len(fine.m_avg))*fine.dt
    try:
        ax.plot(x_axis_vals, fine.m_avg, color='b')
    except ValueError:
        fine.debug(0, "Plot axis value mismatch:")
        fine.debug(0, "N={}".format(len(x_axis_vals)))
        fine.debug(0, "m_avg={}".format(len(fine.m_avg)))
        return
    ax.set_title("Velocity Autocorrelation")
    ax.set_xlabel("Time Lag (s)")
    ax.set_ylabel("$R_{ii}$")
    ax.set_ylim(-1,1)
    # Write to file if write_path given
    if (write_path != None):
        folder = os.path.dirname(write_path)
        if not (os.path.exists(folder)):
            fine.debug(0, "Graph could not be written to file.")
            fine.debug(0, "{} does not exist".format(folder))
        else:
            fig.savefig(write_path)
    if (show):
        plt.show()

def mean_plots2(src, window=80, show=True, write_path=None):
    fine_src = src
    VelocityData.plot_config()
    fig, ax = plt.subplots(1)
    
    fine = VelocityData(fine_src)
    #fine.resize_data(-8000)
    fine.m_avg = pd.rolling_mean(fine.up, window=window)
    weights = np.repeat(1.0, window) / window
    fine.m_avg = np.convolve(fine.up, weights, 'full')
    
    x_axis_vals = np.arange(len(fine.m_avg))*fine.dt
    try:
        ax.plot(x_axis_vals, fine.m_avg, color='b')
    except ValueError:
        fine.debug(0, "Plot axis value mismatch:")
        fine.debug(0, "N={}".format(len(x_axis_vals)))
        fine.debug(0, "m_avg={}".format(len(fine.m_avg)))
        return
    ax.set_title("Velocity Autocorrelation")
    ax.set_xlabel("Time Lag (s)")
    ax.set_ylabel("$R_{ii}$")
    ax.set_ylim(-1,1)
    # Write to file if write_path given
    if (write_path != None):
        folder = os.path.dirname(write_path)
        if not (os.path.exists(folder)):
            fine.debug(0, "Graph could not be written to file.")
            fine.debug(0, "{} does not exist".format(folder))
        else:
            fig.savefig(write_path)
    if (show):
        plt.show()

if __name__ == "__main__":
    main()