
from util import loadAndInterpolateEnsemble, setupMapProjection, decompressVariable, goshen_1km_proj, goshen_1km_gs, probMatchMean, drawPolitical, flux_boxes

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib.ticker import FixedFormatter, FixedLocator
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon

import numpy as np
import Nio as nio

from computeQuantities import computeReflectivity

import glob
from operator import mul
from itertools import izip
import gc
import cPickle
import argparse

def plotConfusion(confusion, map, grid_spacing, title, file_name, inset=None, fudge=16):
    pylab.figure()
    axmain = pylab.axes()
    gs_x, gs_y = grid_spacing

    confusion_cmap_dict = {
        'red':(
            (0.0,  0.5, 0.5), 
            (0.25, 0.5, 1.0),
            (0.5,  1.0, 0.0),
            (0.75, 0.0, 1.0),
            (1.0,  1.0, 1.0),
        ),
        'green':(
            (0.0,  0.5, 0.5), 
            (0.25, 0.5, 0.0),
            (0.5,  0.0, 1.0),
            (0.75, 1.0, 1.0),
            (1.0,  1.0, 1.0),
        ),
        'blue':(
            (0.0,  0.5, 0.5), 
            (0.25, 0.5, 1.0),
            (0.5,  1.0, 0.0),
            (0.75, 0.0, 1.0),
            (1.0,  1.0, 1.0),
        ),
    }
    confusion_cmap = LinearSegmentedColormap('confusion', confusion_cmap_dict, 256)

    nx, ny = confusion.shape
    xs, ys = np.meshgrid( gs_x * np.arange(nx), gs_y * np.arange(ny) )

    pylab.pcolormesh(xs, ys, confusion, cmap=confusion_cmap, vmin=0, vmax=3)

    tick_locs = [ 0.375 + 0.75 * l for l in  range(4) ]
    tick_labels = [ "Correct\nNegative", "False\nAlarm", "Miss", "Hit" ]
    bar = pylab.colorbar()
    bar.locator = FixedLocator(tick_locs)
    bar.formatter = FixedFormatter(tick_labels)
    pylab.setp(pylab.getp(bar.ax, 'ymajorticklabels'), fontsize='large')
    bar.update_ticks()

    drawPolitical(map, scale_len=25)

    if inset:
        lb_y, lb_x = [ b.start for b in inset ]
        ub_y, ub_x = [ b.stop + fudge for b in inset ]

        inset_exp = (slice(lb_y, ub_y), slice(lb_x, ub_x))

        axins = zoomed_inset_axes(pylab.gca(), 2, loc=4)
        pylab.sca(axins)

        pylab.pcolormesh(xs[inset_exp], ys[inset_exp], confusion[inset_exp], cmap=confusion_cmap, vmin=0, vmax=3)
        drawPolitical(map)

        pylab.xlim([lb_x * gs_x, (ub_x - 1) * gs_x])
        pylab.ylim([lb_y * gs_y, (ub_y - 1) * gs_y])

        mark_inset(axmain, axins, loc1=1, loc2=3, fc='none', ec='k')

    pylab.sca(axmain)
    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def ETS(forecast, observation, threshold):
    hits = (forecast >= threshold) & (observation >= threshold)
    misses = (forecast < threshold) & (observation >= threshold)
    false_alarms = (forecast >= threshold) & (observation < threshold)
    correct_negatives = (forecast < threshold) & (observation < threshold)

    num_hits = hits.sum()
    num_misses = misses.sum()
    num_fa = false_alarms.sum()
    num_cn = correct_negatives.sum()

    total = num_hits + num_misses + num_fa + num_cn

    hits_random = float(num_hits + num_misses) * (num_hits + num_fa) / total
    ets = float(num_hits - hits_random) / (num_hits + num_misses + num_fa - hits_random)

    confusion = np.zeros(observation.shape, dtype=np.int32)

    confusion = np.where(correct_negatives, 0, confusion)
    confusion = np.where(false_alarms, 1, confusion)
    confusion = np.where(misses, 2, confusion)
    confusion = np.where(hits, 3, confusion)

    return ets, confusion

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp', dest='exp_name', required=True)
    ap.add_argument('--threshold', dest='threshold', type=int, default=20)

    args = ap.parse_args()

    bounds = (slice(100, 180), slice(90, 170))
    radar_elev, radar_lat, radar_lon = 1883, 41.151944, -104.806111
    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs)
    threshold = args.threshold
    exp_name = args.exp_name
    img_dir = "images-%s/ets_%ddBZ" % (exp_name, threshold)

    map = Basemap(**proj)
    radar_x, radar_y = map(radar_lon, radar_lat)

    obs_base = "hdf/KCYS/1km/goshen.hdfrefl2d"
    obs_times = np.array([ int(f[-6:]) for f in glob.glob("%s*" % obs_base) ])
    fcst_files = glob.glob("/caps1/tsupinie/1km-control-%s/ena???.hdf014[47]00" % exp_name)
    fcst_files.extend(glob.glob("/caps1/tsupinie/1km-control-%s/ena???.hdf01[5678]?00" % exp_name))

    ens_refl, ens_members, ens_times = loadAndInterpolateEnsemble(fcst_files, ['pt', 'p', 'qr', 'qs', 'qh'], computeReflectivity, "/caps1/tsupinie/1km-control-20120712/ena001.hdfgrdbas", 
        {'z_base':radar_elev, 'y_base':radar_y, 'x_base':radar_x, 'elev_angle':0.5}, agl=False, wrap=True)#, aggregator=lambda x: np.mean(x, axis=0))

#   ens_refl, ens_members, ens_times = loadAndInterpolateEnsemble(fcst_files, ['pt', 'p', 'qr', 'qs', 'qh'], computeReflectivity, "/caps1/tsupinie/1km-control-20120712/ena001.hdfgrdbas", 
#       {'z_base':radar_elev, 'y_base':radar_y, 'x_base':radar_x, 'elev_angle':0.5}, agl=False, wrap=True)

#   ens_refl_mean = ens_refl.mean(axis=0)

    refl_ens_mean = probMatchMean(ens_refl)

    bounds_rev = [ slice(None), slice(None) ]
    bounds_rev.extend(bounds[::-1])
    bounds_rev = tuple(bounds_rev)

#   refl_ens_mean = refl_ens_mean[bounds_rev[1:]]
#   ens_refl = ens_refl[bounds_rev]

    all_ets = np.empty((len(ens_members), len(ens_times)), dtype=np.float32)
    all_ets_mean = np.empty((len(ens_times), ), dtype=np.float32)

    all_confusion = np.empty(ens_refl.shape, dtype=np.int32)
    all_confusion_mean = np.empty(refl_ens_mean.shape, dtype=np.int32)

    for wdt, time in enumerate(ens_times):
        idx = np.argmin(np.abs(obs_times - time))
        if obs_times[idx] > time and idx > 0:
            idx -= 1

        bounds_obs = [0]
        bounds_obs.extend(bounds[::-1])
        bounds_obs = tuple(bounds_obs)

        obs_file_name = "%s%06d" % (obs_base, obs_times[idx])
        obs_hdf = nio.open_file(obs_file_name, mode='r', format='hdf')
        obs_refl = obs_hdf.variables['refl2d'][0] #[bounds_obs]

        all_ets_mean[wdt], all_confusion_mean[wdt] = ETS(refl_ens_mean[wdt], obs_refl, threshold)

        gs_x, gs_y = goshen_1km_gs
        for lde, member in enumerate(ens_members):
            all_ets[lde, wdt], all_confusion[lde, wdt] = ETS(ens_refl[lde, wdt], obs_refl, threshold)

#           nx, ny = ens_refl[lde, wdt].shape
#           xs, ys = np.meshgrid( gs_x * np.arange(nx), gs_y * np.arange(ny) )
#           pylab.clf()
#           pylab.contourf(xs, ys, ens_refl[lde, wdt], levels=np.arange(10, 80, 10))
#           pylab.colorbar()
#           pylab.savefig("sweep_interp_%s_%06d.png" % (member, time))

        nx, ny = refl_ens_mean[wdt].shape
        xs, ys = np.meshgrid( gs_x * np.arange(nx), gs_y * np.arange(ny) )
        pylab.clf()
        pylab.contourf(xs, ys, refl_ens_mean[wdt], levels=np.arange(10, 80, 10))
        pylab.colorbar()
        pylab.savefig("%s/sweep_interp_mean_%06d.png" % (img_dir, time))

    cPickle.dump(all_ets_mean, open("%s_%ddBZ.pkl" % (exp_name, threshold), 'w'), -1)

    time_mean_ets = all_ets.mean(axis=1)
    sort_mean_idxs = np.argsort(time_mean_ets)

    pylab.clf()
    for lde, member in enumerate(ens_members):
        print sort_mean_idxs[lde] + 1, time_mean_ets[sort_mean_idxs[lde]]
        pylab.plot(ens_times, all_ets[lde], 'r-', lw=0.75)

    pylab.plot(ens_times, all_ets_mean, 'k-', lw=1.5)

    y_lb, y_ub = pylab.ylim()
    pylab.plot([14400, 14400], [y_lb, y_ub], 'k--', lw=0.5)

    pylab.ylim([y_lb, y_ub])
    pylab.xlim([10800, 18000])

    pylab.xlabel("Time (s)")
    pylab.ylabel("ETS")

    pylab.savefig("%s/ets_swath_mm.png" % img_dir)
    pylab.close()

    for wdt, time in enumerate(ens_times):
        fudge = 16
        if threshold == 20:
            fudge = 32

        plotConfusion(all_confusion_mean[wdt], map, goshen_1km_gs, "Confusion for Reflectivity of the Ensemble Mean at time %06d" % time, "%s/confusion_mean_%06d.png" % (img_dir, time), inset=flux_boxes[exp_name][wdt], fudge=fudge)

 #      for lde, member in enumerate(ens_members):
 #          plotConfusion(all_confusion[lde, wdt], map, goshen_1km_gs, "Confusion for Reflectivity of Member %s at time %06d" % (member, time), "%s/confusion_ena%s_zoom_%06d.png" % (img_dir, member, time))

        gc.collect()

    return

if __name__ == "__main__":
    main()
