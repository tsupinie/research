
import cPickle
from datetime import datetime, timedelta
import glob
from itertools import izip

import numpy as np
import Nio as nio

from util import setupMapProjection, goshen_1km_proj, goshen_1km_gs, loadObs, drawPolitical

import pylab
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon

def plotObservations(obs, map, title, file_name, refl=None):
    pylab.clf()

    if refl is not None:
        nx, ny = refl.shape
        xs, ys = np.meshgrid(np.arange(nx) * 1000, np.arange(ny) * 1000)
        map.contour(xs, ys, refl, levels=np.arange(20, 80, 20), colors='k')

    obs_x, obs_y = map(obs['longitude'], obs['latitude'])

    for ob, x, y in zip(obs, obs_x, obs_y):
        pylab.plot(x, y, 'ko', ms=3.)
        pylab.text(x - 500, y + 500, "%4.1f" % ob['temp'], size='small', va='bottom', ha='right', color='r') 
        pylab.text(x - 500, y - 500, "%4.1f" % ob['dewp'], size='small', va='top', ha='right', color='g') 
        pylab.text(x + 500, y - 500, ob['id'], size='small', va='top', ha='left', color='#808080')

    good_idxs = np.where((obs['wind_spd'] >= 0.) & (obs['wind_dir'] >= 0.))

    u = -obs['wind_spd'] * np.sin(obs['wind_dir'] * np.pi / 180) * 1.94
    v = -obs['wind_spd'] * np.cos(obs['wind_dir'] * np.pi / 180) * 1.94
    u_rot, v_rot = map.rotate_vector(u, v, obs['longitude'], obs['latitude'])
    pylab.barbs(obs_x[good_idxs], obs_y[good_idxs], u_rot[good_idxs], v_rot[good_idxs])

    drawPolitical(map, scale_len=5)

    pylab.title(title)
    pylab.savefig(file_name)
    return

def plotObservationsComposite(obs, map, title, file_name):
    pylab.clf()
    colors = [ 'r', 'g', 'b', 'c', 'm', '#660099', '#ff9900', '#006666' ]

    ob_ids = np.unique1d(obs['id'])
    ttu_label = False

    for ob_id in ob_ids:
        ob_idxs = np.where(obs['id'] == ob_id)[0]
        these_obs = obs[ob_idxs]
        ob_xs, ob_ys = map(these_obs['longitude'], these_obs['latitude'])

        if ob_id[0] == "P":
            ob_num = int(ob_id[1]) - 1
            pylab.plot(ob_xs, ob_ys, 'o', mfc=colors[ob_num], mec=colors[ob_num], ms=3, label="NSSL MM (%s)" % ob_id)
        else:
            if not ttu_label:
                label = "TTU Sticknet"
                ttu_label = True
            else:
                label = None

            pylab.plot(ob_xs[0], ob_ys[0], 'ko', ms=3, label=label)

    drawPolitical(map, scale_len=5)

    pylab.legend(loc=3, numpoints=1, prop={'size':'medium'})
    pylab.title(title)
    pylab.savefig(file_name)
    return

def gatherObservations(obs, target_times):
    thinned_obs = {}
    stn_ids = np.unique(obs['id'])

    for idx, target in enumerate(target_times):
        if idx == 0:
            before_time = target - (target_times[idx + 1] - target) / 2
        else:
            before_time = target - (target - target_times[idx - 1]) / 2

        if idx == len(target_times) - 1:
            after_time = target + (target - target_times[idx - 1]) / 2
        else:
            after_time = target + (target_times[idx + 1] - target) / 2

        keep_obs = np.empty((0,), dtype=obs.dtype)

        for stn in stn_ids:
            stn_idxs = np.where(obs['id'] == stn)[0]

            time_differences = np.abs(obs['time'][stn_idxs] - target)
            closest_time_idx = stn_idxs[np.where(time_differences == np.min(time_differences))[0][0]]
            closest_time = obs['time'][closest_time_idx]


            if before_time < closest_time and closest_time < after_time and np.all((np.abs(keep_obs['latitude'] - obs['latitude'][closest_time_idx]) >= 0.01) | (np.abs(keep_obs['longitude'] - obs['longitude'][closest_time_idx]) >= 0.01)):
                keep_obs.resize(keep_obs.shape[0] + 1)
                keep_obs[-1] = tuple(obs[closest_time_idx])

        thinned_obs[target] = keep_obs
    return thinned_obs

def main():
    _epoch_time = datetime(1970, 1, 1, 0, 0, 0)
    _initial_time = datetime(2009, 6, 5, 18, 0, 0) - _epoch_time
    _initial_time = (_initial_time.microseconds + (_initial_time.seconds + _initial_time.days * 24 * 3600) * 1e6) / 1e6
    _target_times = [ 1800, 3600, 5400, 7200, 9000, 10800, 11100, 11400, 11700, 12000, 12300, 12600, 12900, 13200, 13500, 13800, 14100, 14400,
        14700, 15000, 15300, 15600, 15900, 16200, 16500, 16800, 17100, 17400, 17700, 18000 ]

    inflow_wd_lbound, inflow_wd_ubound = (100, 240)

#   bounds = (0, slice(90, 210), slice(40, 160))
#   bounds = (0, slice(100, 180), slice(90, 170))
    bounds = (0, slice(115, 140), slice(120, 145))
    rev_bounds = [ 0 ]
    rev_bounds.extend(bounds[2:0:-1])
    rev_bounds = tuple(rev_bounds)

    refl_base = "hdf/KCYS/1km/goshen.hdfrefl2d"
    refl_times = np.array([ int(f[-6:]) for f in glob.glob("%s??????" % refl_base) ])
    refl_keep_times = []
    refl_data = {}

    for tt in _target_times:
        idx = np.argmin(np.abs(refl_times - tt))
        if refl_times[idx] > tt and idx > 0:
            idx -= 1

        file_name = "%s%06d" % (refl_base, refl_times[idx])
        hdf = nio.open_file(file_name, mode='r', format='hdf')
        refl_keep_times.append(refl_times[idx])
        refl_data[refl_times[idx]] = hdf.variables['refl2d'][rev_bounds]

    _proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs, bounds=bounds[1:])
#   _proj['resolution'] = 'h' 
    map = Basemap(**_proj)

    ttu_sticknet_obs = cPickle.load(open("ttu_sticknet.pkl", 'r'))
    psu_straka_obs = cPickle.load(open("psu_straka_mesonet.pkl", 'r'))

    all_obs = loadObs(['ttu_sticknet.pkl', 'psu_straka_mesonet.pkl'], [ _epoch_time + timedelta(seconds=(_initial_time + t)) for t in _target_times ],  map, (goshen_1km_proj['width'], goshen_1km_proj['height']), round_time=False)

    partitioned_obs = gatherObservations(all_obs, [ _initial_time + t for t in _target_times ])
    for time, refl_time in zip([ _initial_time + t for t in _target_times], refl_keep_times):
        time_str = (_epoch_time + timedelta(seconds=time)).strftime("%d %B %Y %H%M UTC")

        plot_obs = partitioned_obs[time]  #all_obs[np.where(all_obs['time'] == time)]

        inflow_idxs = np.where((plot_obs['wind_dir'] >= inflow_wd_lbound) & (plot_obs['wind_dir'] <= inflow_wd_ubound))[0]
        outflow_idxs = np.array([ idx for idx in range(plot_obs['id'].shape[0]) if idx not in inflow_idxs ])
#       print "Time = %s" % time_str
#       if len(inflow_idxs) > 0:
#           print "  Inflow stations:", plot_obs[time]['id'][inflow_idxs]
#       else:
#           print "  Inflow stations: [None]"

#       if len(outflow_idxs) > 0:
#           print "  Outflow stations:", plot_obs[time]['id'][outflow_idxs]
#       else:
#           print "  Outflow stations: [None]"

#       print time_str, obs['temp'], obs['dewp'], obs['wind_spd'], obs['wind_dir']

        title = "All MM observations at %s" % time_str
        file_name = "mm_obs_%06d.png" % (time - _initial_time)

        plotObservations(plot_obs, map, title, file_name, refl=refl_data[refl_time])

#   keep_obs = isolateObsTimes(all_obs, [ _epoch_time + timedelta(seconds=(_initial_time + t)) for t in _target_times ])
#   keep_obs.sort(order='time')

#   keep_obs = thinObs(keep_obs, map, goshen_1km_proj['width'], goshen_1km_proj['height'])
    return

if __name__ == "__main__":
    main()
