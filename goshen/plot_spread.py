
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

from dataload import loadEnsemble
from arpsmodelobs import ARPSModelObsFile
from temporal import goshen_3km_temporal, goshen_1km_temporal
from grid import goshen_3km_grid, goshen_1km_grid
from computeQuantities import theta2Temperature

import cPickle

def plotSpread(var, grid, color_bound, title, file_name, refl=None):
    pylab.figure()

    xs, ys = grid.getXY()
    bounds = grid.getBounds()

    pylab.contourf(xs, ys, var[bounds], levels=np.linspace(0, color_bound, 11))
    pylab.colorbar()

    if refl is not None:
        pylab.contour(xs, ys, refl[bounds], colors='k', levels=np.arange(10, 80, 10))

    grid.drawPolitical()

    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def generateBoundary(ens, bounds):
    y_size = bounds[0].stop - bounds[0].start 
    x_size = bounds[1].stop - bounds[1].start

    boundary_spread = {}
    for field in ens.dtype.fields.iterkeys():
        std = ens[field][(slice(None), slice(None)) + bounds].std(axis=0, ddof=1)
        print std.shape

        boundary = np.zeros((std.shape[0], 2 * (y_size + x_size)), dtype=float)
        boundary[:, :x_size]                                 = std[:, 0,  :]
        boundary[:, x_size:(x_size + y_size)]                = std[:, :, -1]
        boundary[:, (x_size + y_size):(2 * x_size + y_size)] = std[:, -1, :]
        boundary[:, (2 * x_size + y_size):]                  = std[:, :,  0]

        boundary_spread[field] = boundary.mean(axis=1)

    cPickle.dump(boundary_spread, open("spread_bnd_%s.pkl" % exp_name, 'w'), -1)
    return

def computeSfc(**kwargs):
    sfc = np.empty(kwargs['pt'].shape, dtype=[('u', np.float32), ('v', np.float32), ('t', np.float32), ('qv', np.float32)])
    sfc['u'] = kwargs['u']
    sfc['v'] = kwargs['v']
    sfc['t'] = theta2Temperature(**kwargs) 
    sfc['qv'] = kwargs['qv']

    return sfc

def main():
    temp = goshen_1km_temporal(start=14400) #start=10800, end=14400)
    grid = goshen_1km_grid() #bounds=(slice(242, 327), slice(199, 284)))
    bounds = grid.getBounds()
    ubounds = {'u':10, 'v':10, 't':3, 'qv':0.0035}    

    base_path = "/caps2/tsupinie/"
    exp_name = "1kmf-mult=1.03"
    exp_path = "%s%s" % (base_path, exp_name)
    n_ens_members = 40
    ens = loadEnsemble(exp_path, n_ens_members, temp.getTimes(), (['u', 'v', 'pt', 'p', 'qv'], computeSfc), points={'sigma':2}, agl=True) #, fcst=True)

    for wdt, (time, time_str) in enumerate(zip(temp, temp.getStrings("%d %B %Y %H%M UTC"))):
        try:
            mo = ARPSModelObsFile("%s/KCYSan%06d" % (exp_path, time))
        except AssertionError:
            mo = ARPSModelObsFile("%s/KCYSan%06d" % (exp_path, time), mpi_config=(2, 12))
 
        for field in ens.dtype.fields.iterkeys():
            std = ens[field][:, wdt, ...].std(axis=0, ddof=1)
            plotSpread(std, grid, ubounds[field], "Spread in %s at %s" % (field, time_str), "spread_%s_%s_%06d.png" % (field, exp_name, time), refl=mo['Z'][0])

    return

if __name__ == "__main__":
   main()
