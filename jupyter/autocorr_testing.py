import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
from velocity_stats import VelocityData

DEBUG = 3

def main():
    N = 800
    dt = 0.00025
    step = 1

    sum_fine = True
    fourier_fine = False
    sum_coarse = False
    fourier_coarse = False
    resize = False

    if sum_fine or fourier_fine:
        fine = VelocityData('./postProcessing/center_probe/290.5/U')
        fine.step = step
        if resize:
            resize_data(fine, -N)

    if sum_coarse or fourier_coarse:
        coarse = VelocityData('./postProcessing_60/center_probe/50/U')
        coarse.step = step
        if resize:
            resize_data(coarse, -N)

    if sum_fine:
        fine.Rxx = autocorr(fine.u, N)
        fine.Ryy = autocorr(fine.v, N)
        fine.Rzz = autocorr(fine.w, N)
        make_plot(fine, N, dt)
    if sum_coarse:
        coarse.Rxx = autocorr(coarse.u, N)
        coarse.Ryy = autocorr(coarse.v, N)
        coarse.Rzz = autocorr(coarse.w, N)
        make_plot(coarse, N, dt)

    if fourier_fine:
        fine.Rxx = autocorr_f(fine.u, N)
        fine.Ryy = autocorr_f(fine.v, N)
        fine.Rzz = autocorr_f(fine.w, N)
        make_plot(fine, N, dt)
    if fourier_coarse:
        coarse.Rxx = autocorr_f(coarse.u, N)
        coarse.Ryy = autocorr_f(coarse.v, N)
        coarse.Rzz = autocorr_f(coarse.w, N)
        make_plot(coarse, N, dt)

    

    #fine.autocorr_sum_plot(lag=8000, step=1, debug=3)
    #coarse.autocorr_sum_plot(lag=8000, step=1, debug=3)

def debug(lvl, msg):
    if not (lvl > DEBUG):
        print(msg)

def reset_data(obj):
    obj.ueff, obj.ufluct = obj.raw, obj.raw - obj.raw.mean()
    obj.u = np.array(obj.ueff[0]['ux'])
    obj.v = np.array(obj.ueff[0]['uy'])
    obj.w = np.array(obj.ueff[0]['uz'])
    obj.up = np.array(obj.ufluct[0]['ux'])
    obj.vp = np.array(obj.ufluct[0]['uy'])
    obj.wp = np.array(obj.ufluct[0]['uz'])
    debug(2, "Data reset complete")
    debug(1, "\nx-component")
    debug(1, "Calculated raw mean  : {}".format(obj.u.mean()))
    debug(1, "Fluctuation raw mean : {}".format(obj.up.mean()))
    debug(1, "\ny-component")
    debug(1, "Calculated raw mean  : {}".format(obj.v.mean()))
    debug(1, "Fluctuation raw mean : {}".format(obj.vp.mean()))
    debug(1, "\nz-component")
    debug(1, "Calculated raw mean  : {}".format(obj.w.mean()))
    debug(1, "Fluctuation raw mean : {}".format(obj.wp.mean()))
    debug(1, "\n")

def resize_data(obj, index):
    obj.u = obj.u[::obj.step][index:]
    obj.v = obj.v[::obj.step][index:]
    obj.w = obj.w[::obj.step][index:]
    obj.up = obj.u - obj.u.mean()
    obj.vp = obj.v - obj.v.mean()
    obj.wp = obj.w - obj.w.mean()
    debug(2, "Data modified for step = {}".format(obj.step))
    debug(2, "Data sliced from index {}".format(index))
    debug(1, "\nx-component")
    debug(1, "Calculated sliced mean  : {}".format(obj.u.mean()))
    debug(1, "Fluctuation sliced mean : {}".format(obj.up.mean()))
    debug(1, "\ny-component")
    debug(1, "Calculated sliced mean  : {}".format(obj.v.mean()))
    debug(1, "Fluctuation sliced mean : {}".format(obj.vp.mean()))
    debug(1, "\nz-component")
    debug(1, "Calculated sliced mean  : {}".format(obj.w.mean()))
    debug(1, "Fluctuation sliced mean : {}".format(obj.wp.mean()))
    debug(1, "\n")

def autocorr(obj, N):
    acovs = []
    try:
        assert not obj.size < N
    except:
        print("Data set {} too small for N value {}".format(obj.size, N))
        sys.exit()
    o = obj[-N:]
    for lag in range(N):
        obj_i = o[:(N-lag)]
        obj_im = pd.Series(o).shift(-lag).dropna()
        acov = np.array(obj_i) * np.array(obj_im)
        acovs.append(acov.sum() / N)
    var = acovs[0]
    acorr = np.array(acovs) / var
    return acorr

def autocorr_f(obj, N):
    o = obj[-N:]
    fvix = np.fft.fft(o, n=2*N)
    #fvix = np.fft.fft(obj)
    acfx = fvix * np.conjugate(fvix)
    acfx = np.fft.ifft(acfx)
    #filter out positive lags
    acfx = acfx[:N]
    #filter out real part
    #acfx = np.real(acfx) / self.N
    #return acfx
    acfx = np.real(acfx)
    return acfx/acfx[0]

def make_plot(obj, N, dt):
    VelocityData.plot_config()
    fig, ax = plt.subplots(1)
    x_axis_vals = np.arange(N)*dt
    try:
        ax.plot(x_axis_vals, obj.Rxx, color='b')
        ax.plot(x_axis_vals, obj.Ryy, color='r')
        ax.plot(x_axis_vals, obj.Rzz, color='g')
    except ValueError:
        debug(0, "Plot axis value mismatch:")
        debug(0, "N={}".format(len(x_axis_vals)))
        debug(0, "Rxx={}".format(len(obj.Rxx)))
        debug(0, "Ryy={}".format(len(obj.Ryy)))
        debug(0, "Rzz={}".format(len(obj.Rzz)))
        return
    ax.set_title("Velocity Autocorrelation")
    ax.set_xlabel("Time Lag (s)")
    ax.set_ylabel("$R_{ii}$")
    ax.set_ylim(-1,1)
    plt.show()


if __name__ == "__main__":
    main()

