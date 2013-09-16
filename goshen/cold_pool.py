
from util import loadObs, interpolate, decompressVariable, drawPolitical, probMatchMean
from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from dataload import loadEnsemble

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import Bbox

import Nio as nio

import numpy as np

from computeQuantities import theta2Temperature, qv2Dewpoint, computeReflectivity

import glob
from datetime import datetime, timedelta
from itertools import izip
import argparse
import cPickle

#def drawPolitical(map):
#    map.drawcoastlines(linewidth=1.5)
#    map.drawcountries(linewidth=1.5)
#    map.drawstates(linewidth=1.0)
#
#    if not hasattr(map, 'counties'):
#        map.readshapefile("countyp020", 'counties', linewidth=0.5)
#    else:
#        for county, data in izip(map.counties_info, map.counties):
#            if county['STATE'] in ['NE', 'WY', 'CO']:
#                pylab.gca().add_patch(Polygon(data, ec='k', fc='none', linewidth=0.5))
#    return

_std_cmap_dict = { 
    'red':(
        (0.0,  0.8,   0.8),
        (0.5,  0.0,   0.0),
        (1.0,  1.0,   1.0), 
    ),
    'green':(
        (0.0,  0.0,   0.0),
        (0.5,  0.0,   0.0),
        (1.0,  0.4,   0.4), 
    ),
    'blue':(
        (0.0,  1.0,   1.0),
        (0.5,  0.0,   0.0),
        (1.0,  0.0,   0.0), 
    ),
}
std_cmap = LinearSegmentedColormap('std_deviation', _std_cmap_dict, 256)

def plotComparison(ens_mean, ens_ob_mean, ens_ob_std, obs, ob_locations, refl, grid, levels, cmap, title, file_name):
    max_std = 5.0

    xs, ys = grid.getXY()
    obs_xs, obs_ys = grid(*ob_locations)

    clip_box = Bbox([[0, 0], [1, 1]])

    pylab.figure()

    pylab.contourf(xs, ys, ens_mean, cmap=cmap, levels=levels)
    pylab.colorbar()

    pylab.contour(xs, ys, refl, colors='k', levels=np.arange(20, 80, 20))

    for ob_x, ob_y, ob, ob_mean, ob_std in zip(obs_xs, obs_ys, obs, ens_ob_mean, ens_ob_std):
        color_bin = np.argmin(np.abs(ob - levels))
        if ob > levels[color_bin]: color_bin += 1
        color_level = float(color_bin) / len(levels)

        ob_z_score = (ob - ob_mean) / ob_std

        print "Ob z-score:", ob_z_score, "Ob:", ob, "Ob mean:", ob_mean, "Ob std:", ob_std

        pylab.plot(ob_x, ob_y, 'ko', markerfacecolor=cmap(color_level), markeredgecolor=std_cmap(ob_z_score / max_std), markersize=4, markeredgewidth=1)
#       pylab.text(ob_x - 1000, ob_y + 1000, "%5.1f" % temp_K, ha='right', va='bottom', size='xx-small', clip_box=clip_box, clip_on=True)

    grid.drawPolitical()

    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def getTempDewpRefl(**kwargs):
    tempdewprefl = np.empty(kwargs['pt'].shape, dtype=[('t', np.float32), ('td', np.float32), ('u', np.float32), ('v', np.float32), ('refl', np.float32)])

    tempdewprefl['t'] = theta2Temperature(**kwargs)
    tempdewprefl['td'] = qv2Dewpoint(**kwargs)
    tempdewprefl['u'] = kwargs['u']
    tempdewprefl['v'] = kwargs['v']
    tempdewprefl['refl'] = computeReflectivity(**kwargs)
    return tempdewprefl

def main():
    np.seterr(all='ignore')
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-path', dest='data_path', default="/caps2/tsupinie/")
    ap.add_argument('--exp-name', dest='exp_name', required=True)

    args = ap.parse_args()

    exp_base = args.data_path
    exp_name = args.exp_name
    n_ensemble_members = 40

#   base_time = datetime(2009, 6, 5, 18, 0, 0)
#   epoch = datetime(1970, 1, 1, 0, 0, 0)
#   base_epoch = (base_time - epoch).total_seconds()
#
#   sec_times = np.arange(14400, 18300, 300)
#   times = [ base_time + timedelta(seconds=int(t)) for t in sec_times ]
    temp = goshen_1km_temporal(start=14400)

    bounds = (slice(100, 180), slice(90, 170))
#   bounds = (slice(None), slice(None))

#   proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs, bounds)
#   map = Basemap(**proj)
    grid = goshen_1km_grid(bounds=bounds)

    obs_file_names = ['psu_straka_mesonet.pkl', 'ttu_sticknet.pkl', 'asos.pkl']
    all_obs = loadObs(obs_file_names, temp.getDatetimes(aslist=True), grid, grid.getWidthHeight())

    obs_x, obs_y = grid(all_obs['longitude'], all_obs['latitude'])
    obs_z = all_obs['elevation']

    grdbas_file = "%s/1kmf-%s/ena001.hdfgrdbas" % (exp_base, exp_name)
    grdbas = nio.open_file(grdbas_file, mode='r', format='hdf')
    y_axis = decompressVariable(grdbas.variables['y'])[bounds[1]]
    x_axis = decompressVariable(grdbas.variables['x'])[bounds[0]]

    y_axis = y_axis - y_axis[0]
    x_axis = x_axis - x_axis[0]

#   fcst_files = glob.glob("%s/1km-control-%s/ena???.hdf014[47]00" % (exp_base, exp_name))
#   fcst_files.extend(glob.glob("%s/1km-control-%s/ena???.hdf01[5678]*" % (exp_base, exp_name)))

#   ens, ens_members, ens_times = loadAndInterpolateEnsemble(fcst_files, ['u', 'v', 'pt', 'p', 'qv', 'qr', 'qs', 'qh'], getTempDewpRefl, grdbas_file, 
#       {'z':10}, agl=True, wrap=True)

    ens = loadEnsemble("/caps2/tsupinie/1kmf-%s/" % exp_name, n_ensemble_members, temp.getTimes(), (['u', 'v', 'pt', 'p', 'qv', 'qr', 'qs', 'qh'], getTempDewpRefl), {'sigma':2}, agl=True, wrap=True) 
#   ens = ens[:, :, 2, :, :]

    ens_slice = [ slice(None), slice(None) ]
    ens_slice.extend(bounds[::-1])

    ens_mean = np.empty(ens.shape[1:], dtype=[('t', np.float32), ('td', np.float32), ('u', np.float32), ('v', np.float32)])
    for var in ens_mean.dtype.fields.iterkeys():
        ens_mean[var] = ens[var].mean(axis=0)

    cPickle.dump(ens_mean, open("cold_pool_1kmf-%s.pkl" % exp_name, 'w'), -1)

    ens = ens[tuple(ens_slice)]

    ens_refl = np.maximum(0, ens['refl'].mean(axis=0)) #probMatchMean(ens['refl'])

    ens_obs = np.empty(ens.shape[:1] + all_obs.shape, dtype=ens_mean.dtype)

    for lde in xrange(ens.shape[0]):
        for var in ens_obs.dtype.fields.iterkeys():
            for ob_idx, (ob_x, ob_y) in enumerate(zip(obs_x, obs_y)):
                wdt = temp.getEpochs(aslist=True).index(int(all_obs['nom_time'][ob_idx]))
                ens_obs[var][lde, ob_idx] = interpolate(ens[var][lde, wdt, np.newaxis], {'y':y_axis, 'x':x_axis}, {'y':ob_y, 'x':ob_y})

#   print ens_obs.shape
    ens_obs_std = np.empty(ens_obs.shape[1:], dtype=ens_obs.dtype)
    ens_obs_mean = np.empty(ens_obs.shape[1:], dtype=ens_obs.dtype)

    for var in ens_obs_std.dtype.fields.iterkeys():
        ens_obs_std[var] = ens_obs[var].std(ddof=1, axis=0)
        ens_obs_mean[var] = ens_obs[var].mean(axis=0)

    cPickle.dump(ens_obs_mean, open("cold_pool_obs_1kmf-%s.pkl" % exp_name, 'w'), -1)

#   print ens_obs_std.shape

#   for wdt, (time_sec, time_epoch) in enumerate(zip(temp, temp.getEpochs())):
#       time_ob_idxs = np.where(all_obs['time'] == time_epoch)[0]
#
#       ob_locations = (all_obs[time_ob_idxs]['longitude'], all_obs[time_ob_idxs]['latitude'])
#
#       temp_K = 5. / 9. * (all_obs[time_ob_idxs]['temp'] - 32) + 273.15
#       dewp_K = 5. / 9. * (all_obs[time_ob_idxs]['dewp'] - 32) + 273.15
#
#       wdir = all_obs[time_ob_idxs]['wind_dir']
#       wspd = all_obs[time_ob_idxs]['wind_spd']
#
#       u = -wspd * np.sin(np.radians(wdir))
#       v = -wspd * np.cos(np.radians(wdir))
#
#       print "Plotting temperature ..."
#       plotComparison(ens_mean['t'][wdt], ens_obs_mean['t'][time_ob_idxs], ens_obs_std['t'][time_ob_idxs], temp_K, ob_locations, ens_refl[wdt], grid, np.arange(289., 298., 1.), matplotlib.cm.get_cmap('Blues_r'),
#           "Ensemble Mean/Obs Comparison at Time %s" % time_sec, "cold_pool_t_%s.png" % time_sec)
#       print "Plotting dewpoint ..."
#       plotComparison(ens_mean['td'][wdt], ens_obs_mean['td'][time_ob_idxs], ens_obs_std['td'][time_ob_idxs], dewp_K, ob_locations, ens_refl[wdt], grid, np.arange(277., 290., 1.), matplotlib.cm.get_cmap('YlGn'),
#           "Ensemble Mean/Obs Comparison at Time %s" % time_sec, "cold_pool_td_%s.png" % time_sec)
#
#       print "Plotting u ..."
#       plotComparison(ens_mean['u'][wdt], ens_obs_mean['u'][time_ob_idxs], ens_obs_std['u'][time_ob_idxs], u, ob_locations, ens_refl[wdt], grid, np.arange(-20., 22., 2.), matplotlib.cm.get_cmap('RdBu_r'),
#           "Ensemble Mean/Obs Comparison at Time %s" % time_sec, "cold_pool_u_%s.png" % time_sec)
#       print "Plotting v ..."
#       plotComparison(ens_mean['v'][wdt], ens_obs_mean['v'][time_ob_idxs], ens_obs_std['v'][time_ob_idxs], v, ob_locations, ens_refl[wdt], grid, np.arange(-20., 22., 2.), matplotlib.cm.get_cmap('RdBu_r'),
#           "Ensemble Mean/Obs Comparison at Time %s" % time_sec, "cold_pool_v_%s.png" % time_sec)
    return

if __name__ == "__main__":
    main()
