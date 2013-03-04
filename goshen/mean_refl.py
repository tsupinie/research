
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

from grid import goshen_3km_grid
from util import loadAndInterpolateEnsemble
from computeQuantities import computeReflectivity

import glob

def plotReflectivity(refl, spread, grid, title, file_name):
    pylab.figure()
    xs, ys = grid.getXY()
    bounds = grid.getBounds()

    pylab.contourf(xs, ys, refl[bounds], levels=np.arange(10, 80, 10))
    pylab.contour(xs, ys, spread[bounds], levels=np.arange(0, 22, 2), cmap=matplotlib.cm.get_cmap('gray_r'))
    pylab.colorbar()

    grid.drawPolitical()

    pylab.suptitle(title)
    pylab.savefig(file_name)
    return

def main():
    grid = goshen_3km_grid(bounds=(slice(242, 327), slice(199, 284)))

    fcst_files = glob.glob("/caps1/tsupinie/3km-control-adapt=0.80/ena???.hdf014400") 
    ens_mean_refl, ens_refl, ens_members, ens_times = loadAndInterpolateEnsemble(fcst_files, ['pt', 'p', 'qr', 'qs', 'qh'], computeReflectivity, "/caps1/tsupinie/3km-control/enf001.hdfgrdbas", 
        {'z':1000}, agl=True, wrap=True, aggregator=lambda x: np.mean(x, axis=0))

    ens_mean_refl = np.maximum(0, np.where(np.isnan(ens_mean_refl), 0, ens_mean_refl))
    ens_refl      = np.maximum(0, np.where(np.isnan(ens_refl),      0, ens_refl))

    ens_spread = ens_refl.std(axis=0, ddof=1)

    for wdt, ens_t in enumerate(ens_times):
        plotReflectivity(ens_mean_refl[wdt], ens_spread[wdt], grid, "Mean and Spread in Reflectivity at %s" % ens_t, "refl_anal_mean_spread_%s-adapt=0.80.png" % ens_t)
    return

if __name__ == "__main__":
    main()
