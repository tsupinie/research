
from legacy import loadAndInterpolateEnsemble, goshen_1km_proj, goshen_1km_gs, setupMapProjection, drawPolitical, inflow_stations
from util import loadObs

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap

import Nio as nio

import numpy as np

import glob
from datetime import datetime, timedelta

def main():
    base_time = datetime(2009, 6, 5, 18, 0, 0)
    epoch = datetime(1970, 1, 1, 0, 0, 0)
    base_epoch = (base_time - epoch).total_seconds()
    times_seconds = range(14700, 18300, 300)
    times = [ base_time + timedelta(seconds=t) for t in times_seconds ]

    bounds = (slice(100, 180), slice(90, 170))
    rev_bounds = [ 0 ]
    rev_bounds.extend(bounds[::-1])
    rev_bounds = tuple(rev_bounds)

    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs, bounds=bounds)
    map = Basemap(**proj)

    obs_file_names = ['psu_straka_mesonet.pkl', 'ttu_sticknet.pkl', 'asos.pkl', 'soundings_clip.pkl']
    all_obs = loadObs(obs_file_names, times, map, (goshen_1km_proj['width'], goshen_1km_proj['height']), sounding_obs=['soundings_clip.pkl'])

    refl_base = "hdf/KCYS/1km/goshen.hdfrefl2d"
    refl_times = np.array([ int(f[-6:]) for f in glob.glob("%s??????" % refl_base) ])
    refl_keep_times = []
    refl_data = {}

    for tt in times_seconds:
        idx = np.argmin(np.abs(refl_times - tt))
        if refl_times[idx] > tt and idx > 0:
            idx -= 1

        file_name = "%s%06d" % (refl_base, refl_times[idx])
        hdf = nio.open_file(file_name, mode='r', format='hdf')
        refl_keep_times.append(refl_times[idx])
        refl_data[tt] = hdf.variables['refl2d'][rev_bounds]

    for time, reg in inflow_stations.iteritems():
        pylab.figure()

        gs_x, gs_y = goshen_1km_gs
        nx, ny = refl_data[time].shape
        xs, ys = np.meshgrid(gs_x * np.arange(nx), gs_y * np.arange(ny))
        pylab.contourf(xs, ys, refl_data[time], levels=np.arange(10, 80, 10))

        for region, stations in reg.iteritems():
            if region != 'sounding':
                for station in stations:
                    idxs = np.where((all_obs['id'] == station) & (all_obs['time'] == base_epoch + time))
                    ob_xs, ob_ys = map(all_obs['longitude'][idxs], all_obs['latitude'][idxs])

                    if region == 'inflow': color='r'
                    elif region == 'outflow': color='b'

                    wdir = all_obs['wind_dir'][idxs]
                    wspd = all_obs['wind_spd'][idxs]
                    u = -wspd * np.sin(wdir * np.pi / 180.) * 1.94
                    v = -wspd * np.cos(wdir * np.pi / 180.) * 1.94

                    pylab.plot(ob_xs, ob_ys, "%so" % color)
                    pylab.barbs(ob_xs, ob_ys, u, v)

        drawPolitical(map, scale_len=10)

        pylab.savefig("inflow_stations_%06d.png" % time)
        pylab.close()

    return

if __name__ == "__main__":
    main()
