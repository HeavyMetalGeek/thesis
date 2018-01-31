import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from itertools import cycle

class MeanFields:
    def __init__(self, 
                 meanPath=None,
                 figPath='postProcessing/figs/',
                 kappa=0.41, 
                 z0=0.00003, 
                 H=1, 
                 u_H=10):
        self.meanPath = meanPath
        self.figPath = figPath
        self.kappa = kappa
        self.z0 = z0
        self.H = H
        self.u_H = u_H
        self.uStar = self.u_H*self.kappa/np.log(self.H/self.z0)
        self.setup()
    
    def setup(self):
        if self.meanPath == None:
            # Allow for meanPath to be unspecified for object creation
            return
        elif not os.path.exists(self.meanPath):
            # Ensure meanPath exists
            raise IOError("Path does not exist! Please change meanPath.")
            return
        # Specify default variables to process
        self.physVars = ['U','uu','vv','ww','uw','vw','uv',
                         'R13','wuu','wvv','www','wuw','wvw','wuv',
                         'R11','nu_SGS','pw','ur13']
        # Create data containers
        self.fileDF = pd.DataFrame(index=self.physVars)
        self.fileDF['filenames'] = ['{0}_mean'.format(i) for i in 
                                    self.physVars]
        self.meanDF = pd.DataFrame()
        self.plot_config()
        

    def time_mean(self, min_time=None, max_time=None):
        # Time limits added so that all data can be combined into a master
        # folder and any time segment can be averaged
        for physVar in self.physVars:
            fname = os.path.join(self.meanPath, 
                                 self.fileDF.at[physVar, 'filenames'])
            df = self.read_mean_files(fname)
            if min_time != None:
                df = df[df.index >= min_time]
            else:
                min_time = df.iloc[0].name
            if max_time != None:
                df = df[df.index <= max_time]
            else:
                max_time = df.iloc[-1].name
            self.meanDF[physVar] = df.mean()
        print("Time mean processed for t={} thru {}".format(min_time, 
                                                            max_time))
        
        self.z = self.meanDF.index.values
        self.dz = np.diff(self.z).mean()
        self.dUdz = np.gradient(self.meanDF['U'], self.dz)
        
        self.U = pd.DataFrame(index=self.meanDF.index)
        self.U['simulation'] = self.meanDF['U']
        self.U['log-law'] = self.uStar/self.kappa*np.log(self.z/self.z0)
        self.dUdz += (-np.gradient(self.U['log-law'], self.dz)
                      + self.uStar / self.kappa / self.z)
        
        self.phi = pd.DataFrame(index=self.meanDF.index)
        self.phi['simulation'] = self.kappa*self.z*self.dUdz/self.uStar
        self.phi['log-law'] = 1
        
        self.vel_var = self.meanDF[['uu', 'vv', 'ww']] / self.uStar**2
        self.vel_var['uuTheo'] = (6 * (1-self.vel_var.index.values)**2 
                                  + self.vel_var.index.values 
                                  * self.vel_var.uu.values[-1])
        self.vel_var['vvTheo'] = (3 * (1-self.vel_var.index.values)**2 
                                  + self.vel_var.index.values 
                                  * self.vel_var.vv.values[-1])
        self.vel_var['wwTheo'] = (1-self.vel_var.index.values)**0.5
        
        self.reynolds = self.meanDF[['uw', 'vw', 'uv']] / self.uStar**2
        
        self.stress = self.meanDF[['uw', 'R13']] / self.uStar**2
        self.stress['total stress'] = self.stress.sum(axis=1)
        self.stress['theoretical'] = self.z - 1
        
        self.tke_flux = self.meanDF[['wuu', 'wvv', 'www']] / self.uStar**3
        
        self.stress_flux = self.meanDF[['wuw', 'wvw', 'wuv']] / self.uStar**3
        
        self.sgs = pd.DataFrame(index=self.meanDF.index)
        self.sgs['kSGS'] = 1.5 * self.meanDF['R11'] / self.uStar**2
        self.sgs['nuSGS'] = self.meanDF['nu_SGS'] / self.uStar / self.H
        self.sgs['-R13'] = -self.meanDF['R13'] / self.uStar**2
        
        self.tke_budget = pd.DataFrame(index=self.meanDF.index)
        self.tke_budget['shear prod'] = -self.meanDF['uw'] * self.dUdz
        self.tke_budget['transport'] = (
                -np.gradient(0.5
                             * self.meanDF[['wuu', 'wvv', 'www']].sum(axis=1)
                             + self.meanDF['pw']
                             + self.meanDF['ur13'],
                             self.dz))
        self.tke_budget['turb. SGS dissip.'] = -self.tke_budget.sum(axis=1)
        self.tke_budget *= self.H / self.uStar**3
        self.tke_budget['total SGS dissip.'] = (
                self.tke_budget['turb. SGS dissip.'] 
                + self.H * self.meanDF['R13'] 
                * self.dUdz / self.uStar**3)
        self.tke_budget['model'] = (1 - 1/self.tke_budget.index) / self.kappa
    
    def read_mean_files(self, fname, zScale=1):
        colNames = ['t', 'dt']
        zcols = pd.np.loadtxt(os.path.join(os.path.dirname(fname),         
                                    'hLevelsCell'))
        zcols = zcols/zScale
        colNames += list(zcols)
        useCols = ['t'] + list(zcols)
        df = pd.read_csv(fname, 
                         delimiter=' ', 
                         comment='#', 
                         names=colNames,
                         usecols=useCols,
                         index_col='t')
        df.columns.name = 'z/H'
        return df

    def adjust_time_mean(self):
        self.uStar = (self.U['simulation'].iloc[-1] 
                      * self.kappa 
                      / np.log(self.H/self.z0))
        self.timeMean()
    
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


    def lines_plot(self, df, **kwargs):

        xlabel = kwargs.get('xlabel', 'NONE')
        linestyles = kwargs.get('linestyles', ['-', '--', '-.', ':'])
        colors = kwargs.get('colors', ['k','r','g','c','m','b'])
        markers = kwargs.get('markers', ['s','o','v','D'])
        
        

        y = df.index
        x = df.values
        
        fig,ax = plt.subplots(1)
        markerCycle = cycle(markers)
        colorCycle = cycle(colors)
        linestyleCycle = cycle(linestyles)
        if x.ndim > 1:
            legends=df.columns
            for i in range(x.shape[1]):
                ax.plot(x[:,i], 
                        y, 
                        marker=next(markerCycle), 
                        color=next(colorCycle),
                        ls=next(linestyleCycle), 
                        label=legends[i], 
                        markeredgecolor='none')
        else:
            ax.plot(x,y,'k-')
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(df.index.name)
        ax.legend(loc='best')
        fig.tight_layout()
        
        return fig,ax


    def plot_val(self, *args, display=True):
        figs = []
        axes = []
        if 'U' in args:
            fig, ax = self.lines_plot(self.U,
                                      xlabel='Mean Streamwise Velocity', 
                                      linestyles=['', '--'],
                                      colors=['k'],
                                      markers=['s', ''])
            ax.set_yscale('log')
            # Limits set to match ASCE paper
            ax.set_xlim([4.875,12])
            figs.append(fig)
            axes.append(ax)
            fname = 'U_mean.png'

        if 'phi' in args:
            fig, ax = self.lines_plot(self.phi,
                                      xlabel='$\Phi$',
                                      linestyles=['', '--'],
                                      colors=['k'],
                                      markers=['s', ''])
            # Limits set to match NIST TN1944
            ax.set_xlim([0.5,2.5])
            figs.append(fig)
            axes.append(ax)
            fname = 'phi_mean.png'

        if 'velocity_variance' in args:
            fig, ax = self.lines_plot(self.vel_var,
                                      xlabel='Velocity Variances',
                                      linestyles=['', '', '', '-', '-', '-'],
                                      colors=['k', 'r', 'g'],
                                      markers=['s', 'o', 'v', '', '', ''])
            # Limits set to match NIST TN1944
            ax.set_xlim([0,8])
            ax.legend(['$\\langle u\'u\'\\rangle$', 
                       '$\\langle v\'v\'\\rangle$', 
                       '$\\langle w\'w\'\\rangle$'], 
                      loc='best')
            figs.append(fig)
            axes.append(ax)
            fname = 'variance_mean.png'

        if 'reynolds_stress' in args:
            fig, ax = self.lines_plot(self.reynolds,
                                      xlabel='Mean Reynolds Stresses$/u_*^2$')
            ax.legend(['$\\langle u\'w\'\\rangle$',
                       '$\\langle v\'w\'\\rangle$',
                       '$\\langle u\'v\'\\rangle$'], 
                      loc='best')
            figs.append(fig)
            axes.append(ax)
            fname = 'reStress_mean.png'

        if 'stress' in args:
            fig, ax = self.lines_plot(self.stress,
                                      xlabel='Mean Stresses$/u_*^2$',
                                      linestyles=['', '', '', '-'],
                                      colors=['k', 'r', 'g', 'g'],
                                      markers = ['s', 'o', 'v', ''])
            # Limits set to match NIST TN1944
            ax.set_xlim([-1,0])
            ax.legend(['$\\langle u\'w\'\\rangle$', 
                       '$\\langle \\tau_{13}\\rangle$',
                       'Total'],
                      loc='best')
            figs.append(fig)
            axes.append(ax)
            fname = 'allStress_mean.png'

        if 'tke_flux' in args:
            fig, ax = self.lines_plot(self.tke_flux,
                                      xlabel=('Flux of Velocity '
                                              'Variances$/u_*^3$'))
            # Limits set to match NIST TN1944
            ax.set_xlim([-0.5,2])
            ax.legend(['$\\langle w\'u\'u\'\\rangle$', 
                       '$\\langle w\'v\'v\'\\rangle$', 
                       '$\\langle w\'w\'w\'\\rangle$'], 
                      loc='best')
            figs.append(fig)
            axes.append(ax)
            fname = 'fluxVariance_mean.png'

        if 'stress_flux' in args:
            fig, ax = self.lines_plot(self.stress_flux,
                                      xlabel=('Flux of Reynolds '
                                              'Stresses$/u_*^3$'))
            # Limits set to match NIST TN1944
            ax.set_xlim([-0.7,0.1])
            ax.legend(['$\\langle w\'u\'w\'\\rangle$', 
                       '$\\langle w\'v\'w\'\\rangle$', 
                       '$\\langle w\'u\'v\'\\rangle$'], 
                      loc='best')
            figs.append(fig)
            axes.append(ax)
            fname = 'fluxStress_mean.png'

        if 'sgs' in args:
            fig, ax = self.lines_plot(self.sgs,
                                      xlabel=('SGS Kinetic Energy, '
                                              'Viscosity \& Stress'))
            ax.set_xscale('log')
            ax.set_yscale('log')
            # Limits set to match NIST TN1944
            ax.set_xlim([10**-5,10**2])
            ax.set_ylim([0.003,1])
            ax.legend(['$\\langle k_{SGS}\\rangle$', 
                       '$\\langle \\nu_{SGS}\\rangle$', 
                       '$-\\langle \\tau_{xz}\\rangle$'], 
                      loc='best')
            figs.append(fig)
            axes.append(ax)
            fname = 'sgs_mean.png'

        if 'tke_budget' in args:
            fig, ax = self.lines_plot(self.tke_budget,
                                      xlabel=('TKE budgets normalized '
                                              'by $u_*^3/H$'))
            ax.set_xlim(-200,200)
            figs.append(fig)
            axes.append(ax)
            fname = 'tke_budgets_mean.png'

        if 'sgs_dissipation' in args:
            fig, ax = self.lines_plot(self.tke_budget[['total SGS dissip.', 
                                                      'model']],
                                      xlabel=('Total SGS dissipation '
                                              'normalized by $u_*^3/H$'),
                                      linestyles=['', '-'],
                                      colors=['g', 'k'],
                                      markers=['s', ''])
            ax.set_xlim(-200,50)
            figs.append(fig)
            axes.append(ax)
            fname = 'sgsDissip_mean.png'

        self.figs = figs
        self.axes = axes
        if display:
            plt.show()
        else:
            pass
            #for f in figs:
            #    f.savefig(self.figPath + fname)

    def combine_data(self,
                     src_dir='./postProcessing/averaging/0/', 
                     dest_dir='./postProcessing/averaging/Combo/', 
                     append=False,
                     min_time=None,
                     max_time=None):
        if not os.path.exists(src_dir):
            raise IOError("Path does not exist! "
                          "Please change source directory.")
            return
        if not os.path.exists(dest_dir):
            os.mkdir(dest_dir)
        fileNames = os.listdir(src_dir)
        if append:
            writeMode = 'a'
            fileNames.remove('hLevelsCell')
        else:
            writeMode = 'w'
        colNames = ['t', 'dt']
        zcols = np.loadtxt(os.path.join(os.path.dirname(src_dir), 
                                        'hLevelsCell'))
        colNames += list(zcols)
        for f in fileNames:
            src_fname = os.path.join(os.path.dirname(src_dir), f)
            df = pd.read_csv(src_fname, 
                                delimiter=' ', 
                                comment='#', 
                                names=colNames, 
                                index_col='t')
            if min_time != None:
                df = df[df.index >= min_time]
            if max_time != None:
                df = df[df.index <= max_time]
            dest_fname = os.path.join(os.path.dirname(dest_dir), f)
            df.to_csv(dest_fname,
                        mode=writeMode,
                        sep=' ',
                        float_format='%.8f',
                        header=False)

def combine_all(src_dir='./postProcessing/averaging/', 
                dest_dir='./postProcessing/averaging/all/'):
    if not os.path.exists(src_dir):
        raise IOError("Path does not exist! "
                      "Please change source directory.")
        return
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)
    exclude_dirs = ['all', 'Combo', 'Combo.bak']
    all_dirs = next(os.walk(src_dir))[1]
    for i in exclude_dirs:
        if i in all_dirs:
            all_dirs.remove(i)
            print("Excluding {}".format(os.path.join(src_dir, i+'/')))
    try:
        all_dirs.sort(key=lambda x: float(x))
    except ValueError:
        print("Not all non-time folders were excluded.")
        return
    dirs = []
    for i in all_dirs:
        dirs.append(os.path.join(src_dir, i+'/'))
    
    # Copy data from hLevelsCell within first folder in src_dir
    hLevelsCell = os.path.join(dirs[0], 'hLevelsCell')
    zcols = np.loadtxt(hLevelsCell)
    colNames = ['t', 'dt'] + list(zcols)
    df = pd.read_csv(hLevelsCell,
                     delimiter=' ',
                     comment='#',
                     names=colNames,
                     index_col='t')
    dest_hLevelsCell = os.path.join(dest_dir, 'hLevelsCell')
    df.to_csv(dest_hLevelsCell,
              mode='w',
              sep=' ',
              float_format='%.3f',
              header=False)

    # Initialized to a None to ensure proper operation 
    min_time = None
    # Initialized to a value which will become min_time
    max_time = -1.0
    # Ensures new files created
    writeMode = 'w'
    for d in dirs:
        # List all files to read and exclude hLevelsCell
        print("\nProcessing data from {}".format(d))
        fileNames = os.listdir(d)
        fileNames.remove('hLevelsCell')
        setTimes = True
        for f in fileNames:
            src_fname = os.path.join(d, f)
            df = pd.read_csv(src_fname, 
                             delimiter=' ', 
                             comment='#', 
                             names=colNames, 
                             index_col='t')
            # Eliminate overlapping time data by truncating data at the
            # beginning of subsequent time files
            if setTimes == True:
                # Set min_time equal to previous directory max_time
                # and max_time equal to current data max time
                min_time, max_time = max_time, df.iloc[-1].name
                setTimes = False
                if min_time > df.iloc[0].name:
                    print("Processing times {} thru {}".format(min_time,
                                                           max_time))
                else:
                    print("Processing times {} thru {}".format(df.iloc[0].name,
                                                               max_time))
            # Truncate current data to times greater 
            # than previous time data max_time
            if not df.iloc[0].name > min_time:
                df = df.loc[df.index > min_time, :]
            
            # Write data to file
            # Write on first time through, append otherwise
            dest_fname = os.path.join(dest_dir, f)
            if writeMode == 'w':
                print("Writing to {}".format(dest_fname), end=' ')
                print("from time {} to {}".format(df.iloc[0].name,
                                                  df.iloc[-1].name))
            if writeMode == 'a':
                print("Appending to {}".format(dest_fname), end=' ')
                print("from {} to {}".format(df.iloc[0].name,
                                             df.iloc[-1].name))
            df.to_csv(dest_fname,
                        mode=writeMode,
                        sep=' ',
                        float_format='%.8f',
                        header=False)
        # Change to append files after initial file set created
        writeMode = 'a'
