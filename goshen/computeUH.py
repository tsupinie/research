
import numpy as np
from scipy.integrate import trapz
from scipy.signal import convolve2d

import matplotlib
matplotlib.use('agg')
import pylab

import argparse
import cPickle

from dataload import loadEnsemble, getAxes
from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from computeQuantities import computeVH
from interp import getInterpFunctions
from util import runConcurrently

def integrateUH(vert_hel, bounds, axes, interp_func):
    updraft_lb, updraft_ub = bounds

    vert_hel_lb = interp_func(vert_hel, axes, {'z':updraft_lb})
    vert_hel_ub = interp_func(vert_hel, axes, {'z':updraft_ub})

    updraft_hel = np.empty(vert_hel_lb.shape, dtype=vert_hel.dtype)

    for idx in np.ndindex(vert_hel_lb.shape):
        full_idx = (slice(None),) + idx

        kdz_lb = np.argmin(np.abs(axes['z'][full_idx] - updraft_lb))
        kdz_ub = np.argmin(np.abs(axes['z'][full_idx] - updraft_ub))

        if axes['z'][(kdz_lb,) + idx] < updraft_lb: kdz_lb += 1
        if axes['z'][(kdz_ub,) + idx] < updraft_ub: kdz_ub += 1

        z_updrft = np.append(np.append(np.array([ updraft_lb ]), axes['z'][(slice(kdz_lb, kdz_ub),) + idx]), updraft_ub)
        hel_updrft = np.append(np.append(np.array([ vert_hel_lb[idx] ]), vert_hel[(slice(kdz_lb, kdz_ub),) + idx]), vert_hel_ub[idx])

        updraft_hel[idx] = trapz(hel_updrft, z_updrft)

    return updraft_hel

def neighborhoodEnsembleProbability(ensemble, threshold):
    kernel = np.array([
        [0., 1., 1., 1., 0.],
        [1., 1., 1., 1., 1.],
        [1., 1., 1., 1., 1.],
        [1., 1., 1., 1., 1.],
        [0., 1., 1., 1., 0.]
    ])

    neighborhood_ensemble = np.empty(ensemble.shape)
    for idx in np.ndindex(ensemble.shape[:-2]):
        neighborhood_ensemble[idx] = convolve2d(ensemble[idx] >= threshold, kernel, mode='same', boundary='symm')

    return neighborhood_ensemble.sum(axis=0) / (neighborhood_ensemble.shape[0] * kernel.sum())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-path', dest='data_path', default="/caps2/tsupinie/")
    ap.add_argument('--exp-name', dest='exp_name', required=True)

    args = ap.parse_args()

    data_path = args.data_path
    exp_name = args.exp_name 
    base_path = "%s/%s/" % (data_path, exp_name)
    n_ens_members = 40
    domain_bounds = (slice(90, 160), slice(90, 160))
    grid = goshen_1km_grid(bounds=domain_bounds)
    temp = goshen_1km_temporal(start=14400)

    tornado_track = zip(*((41.63, -104.383), (41.6134, -104.224)))
    track_xs, track_ys = grid(*[ np.array(list(t)) for t in reversed(tornado_track) ])

    updraft_lb = 0
    updraft_ub = 3000
    updraft_hel_thresh = 125.

    vert_hel = loadEnsemble(base_path, n_ens_members, temp.getTimes(), ([ 'u', 'v', 'w', 'dx', 'dy' ], computeVH), agl=True)

    cPickle.dump(vert_hel[3], open("updraft_hel/VH_%s.pkl" % exp_name, 'w'), -1)

    axes = getAxes(base_path, agl=True)
    interp_func, bounds_func = getInterpFunctions(axes, {'z':updraft_lb})
    
    updraft_hel = np.empty(vert_hel.shape[:2] + vert_hel.shape[-2:], dtype=vert_hel.dtype)

    n_ens_members, n_times = vert_hel.shape[:2]

    for wdt in xrange(n_times):
        print "Integrating time %06d ..." % temp[wdt]

        vh_members = [ vert_hel[(lde, wdt)] for lde in xrange(n_ens_members) ]
        timestep = runConcurrently(integrateUH, vh_members, args=("__placeholder__", (updraft_lb, updraft_ub), axes, interp_func))

        for lde in xrange(n_ens_members):
            updraft_hel[(lde, wdt)] = timestep[lde]

#   print np.nanmax(updraft_hel)
#   argmax_uh = np.unravel_index(np.nanargmax(updraft_hel), updraft_hel.shape)
#   print argmax_uh

    xs, ys = grid.getXY()
    prob_color_map = matplotlib.cm.RdYlBu_r
    prob_color_map.set_under('#ffffff')

    for updraft_hel_thresh in range(75, 425, 25):
        prob_updraft_hel = np.nansum(np.nanmax(updraft_hel, axis=1) >= updraft_hel_thresh, axis=0) / float(n_ens_members)
        nep_updraft_hel = neighborhoodEnsembleProbability(np.nanmax(updraft_hel, axis=1), updraft_hel_thresh)

        pylab.figure()
        pylab.pcolormesh(xs, ys, prob_updraft_hel[domain_bounds], vmin=0.1, vmax=1.0, cmap=prob_color_map)
        pylab.colorbar()

        pylab.plot(track_xs, track_ys, 'mv-', lw=2.5, mfc='k', ms=8)

        grid.drawPolitical()
        pylab.savefig("updraft_hel/UH0-3_%dm2s2_%s.png" % (int(updraft_hel_thresh), exp_name))
        pylab.close()

        pylab.figure()
        pylab.pcolormesh(xs, ys, nep_updraft_hel[domain_bounds], vmin=0.1, vmax=1.0, cmap=prob_color_map)
        pylab.colorbar()

        pylab.plot(track_xs, track_ys, 'mv-', lw=2.5, mfc='k', ms=8)

        grid.drawPolitical()
        pylab.savefig("updraft_hel/NEPcirc_UH0-3_%dm2s2_%s.png" % (int(updraft_hel_thresh), exp_name))
        pylab.close()

        if updraft_hel_thresh == 75:
            cPickle.dump(nep_updraft_hel, open("updraft_hel/NEPcirc_UH0-3_%dm2s2_%s.pkl" % (int(updraft_hel_thresh), exp_name), 'w'), -1)

#   for lde in xrange(n_ens_members):
#       pylab.figure()
#       pylab.pcolormesh(xs, ys, updraft_hel[lde, 0], vmin=25, vmax=125, cmap=prob_color_map)
#       pylab.colorbar()
#       grid.drawPolitical()
#       pylab.savefig("updraft_hel/UH%03d_%s_%06d.png" % (lde + 1, exp_name, 14400))
#       pylab.close()
    return

if __name__ == "__main__":
    main()
