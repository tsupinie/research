
import numpy as np

import Nio as nio

import matplotlib
matplotlib.use('agg')
import pylab

from radarobsfile import RadarObsFile
from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from plot_sounding import plotWinds, plotProfile, plotSkewTBackground
from dataload import loadEnsemble
from computeQuantities import theta2Temperature, qv2Dewpoint

def getSndParams(**kwargs):
    snd_params = np.empty(kwargs['pt'].shape, dtype=[('u', np.float32), ('v', np.float32), ('pt', np.float32), ('p', np.float32), ('qv', np.float32)])
    snd_params['u'] = kwargs['u']
    snd_params['v'] = kwargs['v']
    snd_params['p'] = kwargs['p']
    snd_params['pt'] = kwargs['pt']
    snd_params['qv'] = kwargs['qv']
#   snd_params['t'] = theta2Temperature(**kwargs)
#   snd_params['td'] = qv2Dewpoint(**kwargs)
    return snd_params

def undoAdaptiveInfl(ens_fcst, ens_anal, infl_factor):
    std_fcst = ens_fcst.std(axis=0, ddof=1)
    std_anal = ens_anal.std(axis=0, ddof=1)

    mean_anal = ens_anal.mean(axis=0)
    return (ens_anal - mean_anal) * (std_anal - infl_factor * std_fcst) / (std_anal * (1 - infl_factor)) + mean_anal

def weight(dist, roi, damping_dist):
    if roi > 0:
        z = (dist - roi) / (roi * damping_dist)
    else:
        z = np.zeros(dist.shape, dtype=float)
   
    weights = np.zeros(dist.shape, dtype=float)
    np.where(dist < roi, 1., weights)
    np.where((dist >= roi) & (dist <= ((damping_dist + 1) * roi)), 1. - z, weights)

    return weights

def computeInflWeights(model_z, obs_z, grid, isobs, roi, damping_dist, target=None):
    xs, ys = grid.getXY()
    target_x = xs[target]
    target_y = ys[target]

    inrange = (np.abs(xs - target_x) <= ((damping_dist + 1) * roi)) & (np.abs(ys - target_y) <= ((damping_dist + 1) * roi))
    y_idxs, x_idxs = np.where(inrange)
    x_min, x_max = x_idxs.min(), x_idxs.max() + 1
    y_min, y_max = y_idxs.min(), y_idxs.max() + 1

    horiz_dists = np.zeros(model_z.shape, dtype=float)
    vert_dists = np.zeros(model_z.shape, dtype=float)

    for idx in np.ndindex(model_z[:, y_min:y_max, x_min:x_max].shape):
        grid_idx = tuple(c + o for c, o in zip(idx, [0, y_min, x_min]))

        hdist = np.hypot(xs[grid_idx[1:]] - xs, ys[grid_idx[1:]] - ys)
        hdist = np.where(isobs[:, y_min:y_max, x_min:x_max], hdist, roi * (damping_dist + 2))
        horiz_dists[idx] = hdist.min()
    
        slc = (slice(None), ) + target
        vdist = np.abs(model_z[idx] - obs_z[slc])
        vdist = np.where(isobs[slc], vdist, roi * (damping_dist + 2))
        vert_dists[idx] = vdist.min()

    horiz_weight = weight(horiz_dists, roi, damping_dist)
    vert_weight = weight(vert_dists, roi, damping_dist)
    return horiz_weight * vert_weight

def doMultiplicativeInfl(ensemble, weights, infl_factor):
    for field in ensemble.dtype.fields.iterkeys():
        ens_mean = ensemble[field].mean(axis=0)
        ens_inflt = (ensemble[field] - ens_mean) * infl_factor + ens_mean
        ensemble[field] = np.where(weights[np.newaxis, ...] > 0, weights * ens_inflt, ens_inflt)
    return

def main():
    base_path = "/caps2/tsupinie/1kmf-control/"
    temp = goshen_1km_temporal(start=14400, end=14400)
    grid = goshen_1km_grid()
    n_ens_members = 40

    x_snd, y_snd = grid.getXY(115, 115)
    ens_anal = loadEnsemble(base_path, n_ens_members, temp.getTimes(), (['u', 'v', 'pt', 'p', 'qv'], getSndParams), {'x':x_snd, 'y':y_snd}, fcst=False)
    ens_fcst = loadEnsemble(base_path, n_ens_members, temp.getTimes(), (['u', 'v', 'pt', 'p', 'qv'], getSndParams), {'x':x_snd, 'y':y_snd}, fcst=True)

    robs = RadarObsFile("qc/1km/KCYS.20090605.220000")
    grdbas = nio.open_file("%s/ena001.hdfgrdbas" % base_path, mode='r', format='hdf')
    weights = computeInflWeights(grdbas.variables['zp'], robs.heights, grid, robs['Z'] > 20., 6000, 0)

    ens_mean = np.empty(ens_anal.shape[1:], dtype=ens_anal.dtype)
    ens_preinfl = np.empty(ens_anal.shape, dtype=ens_anal.dtype)

    for field in ens_anal.dtype.fields.iterkeys():
        ens_mean[field] = ens_anal[field].mean(axis=0)
        ens_preinfl[field] = undoAdaptiveInfl(ens_fcst[field], ens_anal[field], 0.9)

    for wdt, (t_ens, time_str) in enumerate(zip(temp, temp.getStrings("%d %B %Y %H%M UTC"))):
        pylab.figure(figsize=(8, 10))
        plotSkewTBackground(pylab.gca())

        for n_ens in xrange(n_ens_members):
            pres_profile = ens_preinfl['p'][n_ens, wdt]
            temp_profile = theta2Temperature(pt=ens_preinfl['pt'][n_ens, wdt], p=ens_preinfl['p'][n_ens, wdt])
            dewp_profile = qv2Dewpoint(qv=ens_preinfl['qv'][n_ens, wdt], p=ens_preinfl['p'][n_ens, wdt])

            if np.any(temp_profile < dewp_profile):
                print "Dewpoint greater than temperature at t=%06d, n=%03d" % (t_ens, n_ens + 1), np.where(temp_profile < dewp_profile)

            plotProfile(temp_profile - 273.15, pres_profile / 100., color='r', linewidth=0.5)
            plotProfile(dewp_profile - 273.15, pres_profile / 100., color='g', linewidth=0.5)

        mean_pres_profile = ens_mean['p'][wdt]
        mean_temp_profile = theta2Temperature(pt=ens_mean['pt'][wdt], p=ens_mean['p'][wdt])
        mean_dewp_profile = qv2Dewpoint(qv=ens_mean['qv'][wdt], p=ens_mean['p'][wdt])

        plotProfile(mean_temp_profile - 273.15, mean_pres_profile / 100., color='k', linewidth=1.5)
        plotProfile(mean_dewp_profile - 273.15, mean_pres_profile / 100., color='k', linewidth=1.5)

        pylab.suptitle("Ensemble Soundings at %s" % time_str)
        pylab.savefig("fcst_snd_1kmf-control_preinfl_%06d.png" % t_ens)
    return

if __name__ == "__main__":
    main()
