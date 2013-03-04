
import Nio as nio
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

import glob
import cPickle

def getTimeFromFilename(filename):
    return int(filename[-6:])

def r_enumerate(L):
    for idx in xrange(len(L) - 1, -1, -1):
        yield idx, L[idx]

def getMaxAndArgmax(refl, times):
    max_refl = refl.max(axis=0)
    argmax_refl = np.argmax(refl, axis=0)
    for wdt, time in r_enumerate(times):
        argmax_refl = np.where(argmax_refl == wdt, time * np.ones(argmax_refl.shape), argmax_refl)

    return max_refl, argmax_refl

def main():
    base_path = "/data6/tsupinie/goshen/hdf/KCYS/1km"
    t_start = 10800
    t_end = 18000

    refl = 0
    filenames = sorted(glob.glob("%s/goshen.hdfrefl2d*" % base_path))
    filenames = [ f for f in filenames if t_start <= getTimeFromFilename(f) and getTimeFromFilename(f) <= t_end ]

    times = [ getTimeFromFilename(f) for f in filenames ]

    for idx, filename in enumerate(filenames):
        hdf = nio.open_file(filename, mode='r', format='hdf')
        if type(refl) == int:
            nz, ny, nx = hdf.variables['refl2d'].shape
            refl = np.empty((len(filenames), ny, nx))

        refl[idx] = hdf.variables['refl2d'][:]
        hdf.close()
    
    cutoff = 14400
    cutoff_idx = np.argmin(np.abs(np.array(times) - cutoff))
    if times[cutoff_idx] < cutoff: cutoff_idx += 1
    max_refl_all,  argmax_refl_all  = getMaxAndArgmax(refl,              times             )
    max_refl_da,   argmax_refl_da   = getMaxAndArgmax(refl[:cutoff_idx], times[:cutoff_idx])
    max_refl_fcst, argmax_refl_fcst = getMaxAndArgmax(refl[cutoff_idx:], times[cutoff_idx:])

    cPickle.dump((max_refl_all,    max_refl_da,    max_refl_fcst),    open("max_obs_refl.pkl",    'w'), -1)
    cPickle.dump((argmax_refl_all, argmax_refl_da, argmax_refl_fcst), open("argmax_obs_refl.pkl", 'w'), -1)

    cmap = matplotlib.cm.jet #RdYlBu_r
    cmap.set_under('#ffffff')
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    pylab.pcolormesh(x, y, max_refl_all, cmap=cmap, vmin=10., vmax=70.)
    pylab.colorbar()

    pylab.savefig("max_refl.png")
    return

if __name__ == "__main__":
    main()
