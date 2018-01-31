
# coding: utf-8

# In[147]:


import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys

get_ipython().run_line_magic('matplotlib', 'inline')


# In[148]:


def plot_config():
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


# ### Process probe data ###

# In[149]:


def process_probes(probe_path):
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


# In[150]:


### Plot Configuration ###
plot_config()

### Processing Parameters ##
# Number of points to process at tail of data #
N = 8000
# Time step
dt = 0.00025
# t0 index
#Nt0 = -N*2-1

### Set which simulation to process ###
# 100 for fine, 60 for coarse
sim = 100
if (sim == 100):
    u = process_probes('./postProcessing/center_probe/290.5/U')
else:
    u = process_probes('./postProcessing_60/center_probe/50/U')

#u = u.iloc[Nt0:]
ueff, up = u, u-u.mean()

u = ueff[0]['ux']
v = ueff[0]['uy']
w = ueff[0]['uz']
fx = up[0]['ux']
fy = up[0]['uy']
fz = up[0]['uz']


# ### Calculating velocity autocorrelation using Fourier transform ###
# $\displaystyle u' = \tilde{u} - U$
# 
# $\displaystyle R_{xx}(m) = \lim_{T\to\infty} \left( \frac{1}{T} \right) \int_0^T f(t) \, f^*(t) \, dt$
# 
# $\displaystyle R_{xx}(m) = \mathcal{F}_x \left[ f_{\nu,i} \; f^*_{\nu,i} \right]$

# In[151]:


Nx = len(fx)
#fvix = np.fft.fft(fx, n=2*Nx)
fvix = np.fft.fft(fx)
acfx = fvix * np.conjugate(fvix)
acfx = np.fft.ifft(acfx)
#filter out positive lags
acfx = acfx[:Nx]
#filter out real part
acfx = np.real(acfx)/Nx

Ny = len(fy)
#fviy = np.fft.fft(fy, n=2*Ny)
fviy = np.fft.fft(fy)
acfy = fviy * np.conjugate(fviy)
acfy = np.fft.ifft(acfy)
#filter out positive lags
acfy = acfy[:Ny]
#filter out real part
acfy = np.real(acfy)/Ny

Nz = len(fz)
#fviz = np.fft.fft(fz, n=2*Nz)
fviz = np.fft.fft(fz)
acfz = fviz * np.conjugate(fviz)
acfz = np.fft.ifft(acfz)
#filter out positive lags
acfz = acfz[:Nz]
#filter out real part
acfz = np.real(acfz)/Nz


# ### Calculating velocity autocorrelation using discrete values ###
# 
# $\displaystyle \hat{R}_{xx}(m) = \left( \frac{1}{N} \right) \sum_{i=0}^{N-m-1} f_i \, f^*_{i+m}$

# In[152]:


Rxx = []
Ryy = []
Rzz = []

#lag = list(range(0, N))
fxi = np.array(fx.iloc[:N])
fyi = np.array(fy.iloc[:N])
fzi = np.array(fz.iloc[:N])

def autocorr(obj, N, lag):
    obji = obj.iloc[-(N+lag):].iloc[:N]
    objim = obj.iloc[-N:]
    acov = np.array(obji) * np.array(objim)
    var = np.array(obji) * np.array(obji)
    return acov.sum() / var.sum()

m = N
for i in range(m):
    Rxx.append(autocorr(fx, N, i))
    Ryy.append(autocorr(fy, N, i))
    Rzz.append(autocorr(fz, N, i))

'''
for m in lag:
    fxim = np.array(fx.iloc[Nt0+m:].iloc[:N])
    Rx = fxi * fxim
    fxvar = (fxi**2).sum()/N
    Rxx.append(Rx.sum()/N/fxvar)
    
    fyim = np.array(fy.iloc[Nt0+m:].iloc[:N])
    #Ry = fyi * fyim.conj()
    #Ryy.append(Ry.sum()/N)
    Ry = fyi * fyim
    fyvar = (fyi**2).sum()/N
    Ryy.append(Ry.sum()/N/fyvar)
    
    fzim = np.array(fz.iloc[Nt0+m:].iloc[:N])
    #Rz = fzi * fzim.conj()
    #Rzz.append(Rz.sum()/N)
    Rz = fzi * fzim
    fzvar = (fzi**2).sum()/N
    Rzz.append(Rz.sum()/N/fzvar)
'''


# ### PLOT VELOCITY FLUCTUATIONS ###

# In[153]:


fig, ax = plt.subplots(1)
print(Nx)
ax.plot(np.arange(Nx)*dt, u, color='b')
ax.axhline(u.mean(), 0, Nx, ls='-.', color='b')
ax.plot(np.arange(Ny)*dt, v, color='r')
ax.axhline(v.mean(), 0, Nx, ls='-.', color='r')
ax.plot(np.arange(Nz)*dt, w, color='g')
ax.axhline(w.mean(), 0, Nx, ls='-.', color='g')
ax.set_title("Velocity Fluctuations")
ax.set_xlabel("Time (s)")
ax.set_ylabel("$u'_i$")
ax.set_xlim(0, np.floor(Nx*dt*100)/100)
if (sim == 100):
    fig.savefig("./figs/velocity_fluctuations.png")
else:
    fig.savefig("./figs/velocity_fluctuations60.png")
plt.show()


# ### PLOT VELOCITY AUTOCORRELATION ###

# In[154]:


fig, ax = plt.subplots(1)
ax.plot(np.arange(0, N)*dt, acfx[-N:], color='b')
ax.plot(np.arange(0, N)*dt, acfy[-N:], color='r')
ax.plot(np.arange(0, N)*dt, acfz[-N:], color='g')
ax.set_title("Velocity Autocorrelation")
ax.set_xlabel("Time Lag (s)")
ax.set_ylabel("$R_{ii}$")
#ax.set_ylim(-1,1)
if (sim == 100):
    fig.savefig('./figs/autocorrelation_fourier.png')
else:
    fig.savefig('./figs/autocorrelation_fourier60.png')
plt.show()


fig, ax = plt.subplots(1)
ax.plot(np.arange(0, N)*dt, Rxx, color='b')
ax.plot(np.arange(0, N)*dt, Ryy, color='r')
ax.plot(np.arange(0, N)*dt, Rzz, color='g')
ax.set_title("Velocity Autocorrelation")
ax.set_xlabel("Time Lag (s)")
ax.set_ylabel("$R_{ii}$")
#ax.set_ylim(-1,1)
if (sim == 100):
    fig.savefig('./figs/autocorrelation_sum.png')
else:
    fig.savefig('./figs/autocorrelation_sum60.png')
plt.show()

from scipy.stats import norm
fig, ax = plt.subplots(1)
ax.plot(u, norm.pdf(u), 'r-', lw=5, alpha=0.6, label='norm pdf')
ax.hist(u, normed=True, histtype='stepfilled', bins='auto')
plt.show()

fig, ax = plt.subplots(1)
ax.plot(v, norm.pdf(v), 'r-', lw=5, alpha=0.6, label='norm pdf')
ax.hist(v, normed=True, histtype='stepfilled', bins='auto')

ax.plot(w, norm.pdf(w), 'r-', lw=5, alpha=0.6, label='norm pdf')
ax.hist(w, normed=True, histtype='stepfilled', bins='auto')
plt.show()

