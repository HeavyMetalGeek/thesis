import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys

class VelocityData:
    def __init__(self, probe_path):
        # Default values
        self.lag = None
        self.N = None
        self.step = 1
        self.dt = 0.00025
        self.probe_dt = 0.00025
        self.debug_lvl = 0
        
        self.path = probe_path
        self.raw = VelocityData.process_probes(self.path)
        self.reset_data()
        
    @classmethod
    def process_probes(cls, probe_path):
        os.system('sed -i_orig \'s/[(,)]//g\' {}'.format(probe_path))
        u = pd.read_csv(probe_path, 
                        delim_whitespace=True, 
                        comment='#', 
                        header=None, 
                        index_col=0)
        num_pts = int(u.shape[1]/3)
        cols = pd.MultiIndex.from_tuples([(x,y) for x in range(num_pts) 
                                          for y in ['ux', 'uy', 'uz']])

        u.columns = cols
        u.index.name = 'Time'
        u.columns.names = ['Probes', 'Component']
        return u

    @classmethod
    def plot_config(cls):
        plt.rc('xtick', labelsize=14)
        plt.rc('ytick', labelsize=14)
        plt.rc('xtick.major', size=8.0, width=1.0)
        plt.rc('xtick.minor', size=4.0, width=1.0)
        plt.rc('ytick.major', size=8.0, width=1.0)
        plt.rc('ytick.minor', size=4.0, width=1.0)
        plt.rc('axes', labelsize=16)
        plt.rc('legend', fontsize=14)
        plt.rc('axes', titlesize=16)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.rc('font', serif='cm')
    
    def debug(self, lvl, msg):
        if not (lvl > self.debug_lvl):
            print(msg)
    
    def process_data(self, probe_path=None, **kwargs):
        if (probe_path != None):
            self.__init__(probe_path)
        self.kwarg_handler(**kwargs)
                
    def reset_data(self):
        self.ueff, self.ufluct = self.raw, self.raw - self.raw.mean()
        self.u = np.array(self.ueff[0]['ux'])
        self.v = np.array(self.ueff[0]['uy'])
        self.w = np.array(self.ueff[0]['uz'])
        self.up = np.array(self.ufluct[0]['ux'])
        self.vp = np.array(self.ufluct[0]['uy'])
        self.wp = np.array(self.ufluct[0]['uz'])
        self.debug(2, "Data reset complete")
        self.debug(1, "Calculated raw mean  : {}".format(self.u.mean()))
        self.debug(1, "Fluctuation raw mean : {}".format(self.up.mean()))

    def resize_data(self, index):
        self.u = self.u[::self.step][index:]
        self.v = self.v[::self.step][index:]
        self.w = self.w[::self.step][index:]
        self.up = self.u - self.u.mean()
        self.vp = self.v - self.v.mean()
        self.wp = self.w - self.w.mean()
        self.debug(2, "Data modified for step = {}".format(self.step))
        self.debug(2, "Data sliced from index {}".format(index))
        self.debug(1, "Calculated sliced mean  : {}".format(self.u.mean()))
        self.debug(1, "Fluctuation sliced mean : {}".format(self.up.mean()))

    def kwarg_handler(self, **kwargs):
        for key, val in kwargs.items():
            if (key == 'N'):
                self.N = val
            if (key == 'lag'):
                self.lag = val
            if (key == 'step'):
                self.step = val
                self.dt = self.step * self.probe_dt
            if (key == 'debug'):
                self.debug_lvl = val
        
    def autocorr(self, summation=True, **kwargs):
        # Process kwargs
        self.kwarg_handler(**kwargs)

        # Start with complete data set
        self.reset_data()

        # Handle missing values for functions
        # Add capability for custom functions
        self.debug(2, "N set to {}".format(self.N))
        self.debug(2, "lag set to {}".format(self.lag))
        self.debug(2, "step set to {}".format(self.step))
        if (summation):
            self.debug(2, "Summation algorithm start")
            try:
                assert self.lag != None
                if (self.N == None):
                    self.N = self.lag
                self.autocorr_sum(self.lag)
                self.debug(2, "Summation computation complete")
            except:
                self.debug(0, "EXCEPTION - parameter undefined: lag")
        else:
            self.debug(2, "Fourier algorithm start")
            try:
                assert self.N != None
                self.autocorr_fourier(self.N)
                self.debug(2, "Fourier computation complete")
            except:
                self.debug(0, "EXCEPTION - parameter undefined: N")

    def autocorr_sum(self, lag):
        if (self.N == None):
            self.N = self.lag
            self.debug(2, "N set equal to lag")
        if (self.lag > self.N):
            self.lag = self.N
            self.debug(2, "Value for lag reduced to {}".format(self.lag))
        self.resize_data(-self.lag)
        self.debug(2, "Data reduced to {} points".format(self.lag))
        
        self.Rxx = self.autocorr_sum_func(self.up)
        self.debug(2, "Rxx sum computation complete")
        self.Ryy = self.autocorr_sum_func(self.vp)
        self.debug(2, "Ryy sum computation complete")
        self.Rzz = self.autocorr_sum_func(self.wp)
        self.debug(2, "Rzz sum computation complete")

    def autocorr_sum_func(self, obj):
        acovs = []
        self.debug(2, "Starting summation loop")
        try:
            for lag in range(self.N):
                obj_i = obj[:self.N-lag]
                obj_im = pd.Series(obj).shift(-lag).dropna()
                acov = np.array(obj_i) * np.array(obj_im)
                acovs.append(acov.sum() / self.N)
            var = acovs[0]
            acorr = np.array(acovs) / var
        except:
            self.debug(0, "EXCEPTION: error within summation loop")
        return acorr

    def autocorr_fourier(self, N):
        self.resize_data(-(self.N))
        
        self.Rxx = self.autocorr_fourier_func(self.up)
        self.debug(2, "Rxx fourier computation complete")
        self.Ryy = self.autocorr_fourier_func(self.vp)
        self.debug(2, "Ryy fourier computation complete")
        self.Rzz = self.autocorr_fourier_func(self.wp)
        self.debug(2, "Rzz fourier computation complete")

    def autocorr_fourier_func(self, obj):
        self.debug(2, "Starting Fourier computation")
        fvix = np.fft.fft(obj, n=2*self.N)
        #fvix = np.fft.fft(obj)
        acfx = fvix * np.conjugate(fvix)
        acfx = np.fft.ifft(acfx)
        #filter out positive lags
        acfx = acfx[:self.N]
        #filter out real part
        #acfx = np.real(acfx) / self.N
        #return acfx
        acfx = np.real(acfx)
        return acfx/acfx[0]

    def velocity_fluctuation_plot(self, show=True, write_path=None, **kwargs):
        self.reset_data()
        self.kwarg_handler(**kwargs)
        if (self.lag == None):
            self.lag = 0
            self.debug(2, "Value for lag set equal to 0")
        if (self.N == None):
            self.N = self.lag
            self.debug(2, "N set equal to lag")
        if (self.lag > self.N):
            self.lag = self.N
            self.debug(2, "Value for lag reduced to {}".format(self.lag))
        self.resize_data(-(self.N + self.lag))
        
        fig, ax = plt.subplots(1)
        try:
            ax.plot(np.arange(self.N + self.lag) * self.dt, self.u, color='b')
            ax.plot(np.arange(self.N + self.lag) * self.dt, self.v, color='r')
            ax.plot(np.arange(self.N + self.lag) * self.dt, self.w, color='g')
            ax.axhline(self.u.mean(), 0, self.N + self.lag, ls='-.', color='k')
            ax.axhline(self.v.mean(), 0, self.N + self.lag, ls='-.', color='k')
            ax.axhline(self.w.mean(), 0, self.N + self.lag, ls='-.', color='k')
        except ValueError:
            self.debug(0, "Plot axis value mismatch:")
            self.debug(0, "N={}".format(len(self.N)))
            self.debug(0, "u={}".format(len(self.u)))
            return
        ax.set_title("Velocity Fluctuations")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("$u_i$")
        # Make sure last value on x-axis is a whole number
        x_lim_max = np.floor(self.N*self.dt*100)/100
        if not (x_lim_max < 1):
            x_lim_max = np.floor(x_lim_max)
        # Make sure last value on x-axis is at graph border
        if (x_lim_max < ax.get_xticks()[-1]):
            ax.set_xticks(ax.get_xticks()[:-1])
            x_lim_max = ax.get_xticks()[-1]
        # Set x-axis limits
        ax.set_xlim(0, x_lim_max)
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

    def autocorr_sum_plot(self, show=True, write_path=None, **kwargs):
        self.autocorr(**kwargs)
        VelocityData.plot_config()
        fig, ax = plt.subplots(1)
        x_axis_vals = np.arange(self.N)*self.dt
        try:
            ax.plot(x_axis_vals, self.Rxx, color='b')
            ax.plot(x_axis_vals, self.Ryy, color='r')
            ax.plot(x_axis_vals, self.Rzz, color='g')
        except ValueError:
            self.debug(0, "Plot axis value mismatch:")
            self.debug(0, "N={}".format(len(x_axis_vals)))
            self.debug(0, "Rxx={}".format(len(self.Rxx)))
            self.debug(0, "Ryy={}".format(len(self.Ryy)))
            self.debug(0, "Rzz={}".format(len(self.Rzz)))
            return
        ax.set_title("Velocity Autocorrelation")
        ax.set_xlabel("Time Lag (s)")
        ax.set_ylabel("$R_{ii}$")
        ax.set_ylim(-1,1)
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

    def autocorr_fourier_plot(self, show=True, write_path=None, **kwargs):
        self.autocorr(summation=False, **kwargs)
        fig, ax = plt.subplots(1)
        x_axis_vals = np.arange(self.N)*self.dt
        try:
            ax.plot(x_axis_vals, self.Rxx, color='b')
            ax.plot(x_axis_vals, self.Ryy, color='r')
            ax.plot(x_axis_vals, self.Rzz, color='g')
        except ValueError:
            self.debug(0, "Plot axis value mismatch:")
            self.debug(0, "N={}".format(len(x_axis_vals)))
            self.debug(0, "Rxx={}".format(len(self.Rxx)))
            self.debug(0, "Ryy={}".format(len(self.Ryy)))
            self.debug(0, "Rzz={}".format(len(self.Rzz)))
            return
        ax.set_title("Velocity Autocorrelation")
        ax.set_xlabel("Time Lag (s)")
        ax.set_ylabel("$R_{ii}$")
        ax.set_ylim(-1,1)
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