
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab
from matplotlib.colors import LinearSegmentedColormap

import cPickle
import glob
from datetime import datetime, timedelta
from math import ceil
import sys

from radarobsfile import RadarObsFile
from grid import goshen_3km_grid

def innovation(obs, forecasts, dindex):
    innovation = (obs - forecasts)[dindex]
    rmsi = np.sqrt((innovation ** 2).mean())
    return rmsi, innovation.mean()

def spread(forecasts, fcst_mean, dindex):
    return np.sqrt(np.mean(np.sum((forecasts - fcst_mean) ** 2, axis=0)[dindex] / (forecasts.shape[0] - 1)))

def matchObsAndForecast(obs_times, forecast_times):
    epoch = datetime(1970, 1, 1, 0, 0, 0)
    obs_round = [ epoch + timedelta(seconds=ceil((o - epoch).total_seconds() / 300.) * 300) for o in obs_times ]

    matches = {}
    for idx, time in enumerate(obs_round):
        if time in forecast_times:
            matches[time] = (idx, np.where(np.array(forecast_times) == time)[0][0])

    return sorted(matches.values(), key=lambda x: x[0])

def main():
    n_ens_members = 40
    radars = ['KCYS', 'KFTG', 'KRIW']

    missing_thd = -90
    stat_thd = 15

    base_time = datetime(2009, 6, 5, 18, 0, 0)
    obs_files = sorted(glob.glob("qc/manual/3km/KCYS.20090605.??????"))
    fcst_files = [ sorted(glob.glob("3km-fixed-radar/eno%03d.pkl??????" % (n_ens + 1))) for n_ens in range(n_ens_members) ]
    fcst_mean_files = sorted(glob.glob("3km-fixed-radar/eomean.pkl??????"))

    obs_times = [ datetime.strptime(of[-15:], "%Y%m%d.%H%M%S") for of in obs_files ]
    fcst_times = [ base_time + timedelta(seconds=int(ff[-6:])) for ff in fcst_files[0] ]

    matches = matchObsAndForecast(obs_times, fcst_times)

    spread_curve_refl = np.empty((len(radars), len(fcst_times)), dtype=float)
    spread_curve_vel = np.empty((len(radars), len(fcst_times)), dtype=float)
    rmsi_curve_refl = np.empty((len(radars), len(fcst_times)), dtype=float)
    rmsi_curve_vel = np.empty((len(radars), len(fcst_times)), dtype=float)
    minnov_curve_refl = np.empty((len(radars), len(fcst_times)), dtype=float)
    minnov_curve_vel = np.empty((len(radars), len(fcst_times)), dtype=float)

    grid = goshen_3km_grid(bounds=(slice(242, 327), slice(199, 284)))

    for obs_idx, fcst_idx in matches:
        print "Loading %s ..." % obs_files[obs_idx]
        obs_data = RadarObsFile(obs_files[obs_idx])
        fcst_data = []

        fcst_mean = cPickle.load(open(fcst_mean_files[fcst_idx], 'r'))

        sys.stdout.write("From %s, loading member " % fcst_files[0][fcst_idx].replace("eno001", "eno???"))
        for n_ens in range(n_ens_members):
            sys.stdout.write("%d ... " % (n_ens + 1))
            sys.stdout.flush()

            fcst_data.append(cPickle.load(open(fcst_files[n_ens][fcst_idx], 'r')))

        sys.stdout.write("\n")
        sys.stdout.flush()

        for rad_idx, radar in enumerate(radars):
            refl_fcst = np.array([ fcst['Z'][rad_idx] for fcst in fcst_data ])
            refl_fcst_mean = fcst_mean['Z'][rad_idx]

            vel_fcst = np.array([ fcst['vr'][rad_idx] for fcst in fcst_data ])
            vel_fcst_mean = fcst_mean['vr'][rad_idx]

            refl_obs = np.transpose(obs_data.reflectivity, [2, 0, 1])
            vel_obs = np.transpose(obs_data.radial_velocity, [2, 0, 1])

            refl_fcst = np.where(refl_fcst <= missing_thd, np.nan, refl_fcst)
            vel_fcst = np.where(vel_fcst <= missing_thd, np.nan, refl_fcst)

            mask = ((refl_obs > stat_thd) | (refl_fcst_mean > stat_thd)) & (refl_obs > missing_thd) & (refl_fcst.mean(axis=0) > missing_thd)
            dindex = np.where(mask)

            spd = spread(refl_fcst, refl_fcst_mean, dindex)
            rmsi, minnov = innovation(refl_obs, refl_fcst_mean, dindex)

            spread_curve_refl[rad_idx, fcst_idx] = spd
            rmsi_curve_refl[rad_idx, fcst_idx] = rmsi
            minnov_curve_refl[rad_idx, fcst_idx] = minnov

            mask = ((vel_obs > stat_thd) | (vel_fcst_mean > stat_thd)) & (vel_obs > missing_thd) & (vel_fcst.mean(axis=0) > missing_thd)
            dindex = np.where(mask)

            spd = spread(vel_fcst, vel_fcst_mean, dindex)
            rmsi, minnov = innovation(vel_obs, vel_fcst_mean, dindex)

            spread_curve_vel[rad_idx, fcst_idx] = spd
            rmsi_curve_vel[rad_idx, fcst_idx] = rmsi
            minnov_curve_vel[rad_idx, fcst_idx] = minnov

        print spread_curve_refl[:, -1], rmsi_curve_refl[:, -1], minnov_curve_refl[:, -1], spread_curve_vel[:, -1], rmsi_curve_vel[:, -1], minnov_curve_vel[:, -1]

#       for kdz in range(14):
#           pylab.clf()
#           xs, ys = grid.getXY()
#           bounds = grid.getBounds()
#           pylab.pcolormesh(xs, ys, np.where(mask[kdz][bounds], refl_obs[kdz][bounds] - refl_fcst_mean[kdz][bounds], np.nan), cmap=matplotlib.cm.get_cmap('RdYlBu'), vmin=-50, vmax=50) 
##          pylab.pcolormesh(xs, ys, refl_fcst_mean[kdz][bounds], cmap=matplotlib.cm.get_cmap('jet'), vmin=10, vmax=80) 
#           pylab.colorbar()
##          pylab.contour(xs + 1500, ys + 1500, refl_obs[kdz][bounds], levels=[ 15. ], colors='k')
##          pylab.pcolormesh(xs, ys, mask[kdz][bounds], cmap=matplotlib.cm.get_map('binary_r'), vmin=0, vmax=1, alpha=0.5)
#           grid.drawPolitical()
#           pylab.savefig("innov_%s_tilt%02d.png" % (fcst_files[0][fcst_idx][-6:], kdz))

    cPickle.dump((spread_curve_refl, spread_curve_vel, rmsi_curve_refl, rmsi_curve_vel, minnov_curve_refl, minnov_curve_vel), open("fcst_stat.pkl", 'w'), -1)

    return

if __name__ == "__main__":
    main()
