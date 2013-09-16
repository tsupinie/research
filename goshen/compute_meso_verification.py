
import numpy as np
from scipy.interpolate import griddata

import matplotlib
matplotlib.use('agg')
import pylab

from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from dataload import loadEnsemble
from util import loadObs
from computeQuantities import theta2Temperature, qv2Dewpoint
from arpsmodelobs import ARPSModelObsFile

import argparse
import cPickle

def getTempDewpRefl(**kwargs):
    tempdewprefl = np.empty(kwargs['pt'].shape, dtype=[('t', np.float32), ('td', np.float32), ('u', np.float32), ('v', np.float32)])

    tempdewprefl['t'] = theta2Temperature(**kwargs)
    tempdewprefl['td'] = qv2Dewpoint(**kwargs)
    tempdewprefl['u'] = kwargs['u']
    tempdewprefl['v'] = kwargs['v']
    return tempdewprefl

def main():
    base_path = "/caps2/tsupinie/"
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp-name', dest='exp_name', required=True)

    args = ap.parse_args()

    n_ens_members = 40
    exp_name = args.exp_name

    bounds_obs = (slice(100, 180), slice(90, 170))
    grid_obs = goshen_1km_grid(bounds=bounds_obs)

    bounds = (slice(None), slice(None))
    grid = goshen_1km_grid(bounds=bounds)

    temp = goshen_1km_temporal(start=14400)

    obs_file_names = ['psu_straka_mesonet.pkl', 'ttu_sticknet.pkl', 'asos.pkl']
    all_obs = loadObs(obs_file_names, temp.getDatetimes(aslist=True), grid_obs, grid_obs.getWidthHeight())
    obs_xy = np.vstack(grid(all_obs['longitude'], all_obs['latitude'])).T

    ens = loadEnsemble("/caps2/tsupinie/%s/" % exp_name, n_ens_members, temp.getTimes(), (['u', 'v', 'pt', 'p', 'qv'], getTempDewpRefl), {'sigma':2}, agl=True, wrap=True) 

    grid_xs, grid_ys = grid.getXY()
    obs_t_verif = []
    for wdt, (time_sec, time_epoch) in enumerate(zip(temp, temp.getEpochs())):

        try:
            mo = ARPSModelObsFile("%s/%s/KCYSan%06d" % (base_path, exp_name, time_sec))
        except AssertionError:
            mo = ARPSModelObsFile("%s/%s/KCYSan%06d" % (base_path, exp_name, time_sec), mpi_config=(2, 12))
        except:
            print "Can't load reflectivity ..."
            mo = {'Z':np.zeros((1, 255, 255), dtype=np.float32)}

        time_ob_idxs = np.where(all_obs['nom_time'] == time_epoch)[0]

        time_obs = all_obs[time_ob_idxs]
        time_obs_xy = obs_xy[time_ob_idxs]

        obs_intrp = griddata(time_obs_xy, 5. / 9. * (time_obs['temp'] - 32) + 273.15, (grid_xs, grid_ys))
        print np.isfinite(obs_intrp).sum()

        pylab.figure()

        pylab.contourf(grid_xs, grid_ys, ens['t'][:, wdt].mean(axis=0)[bounds] - obs_intrp, levels=np.arange(-6, 6.5, 0.5), cmap=matplotlib.cm.get_cmap("RdBu_r"))
        pylab.colorbar()

        pylab.contour(grid_xs, grid_ys, mo['Z'][0][tuple(reversed(bounds))], levels=np.arange(10, 80, 10), colors='k')

        grid.drawPolitical()

        pylab.savefig("obs_verif/obs_%s_t_grid_%06d.png" % (exp_name[5:], time_sec))
        pylab.close()
        obs_t_verif.append(ens['t'][:, wdt].mean(axis=0) - obs_intrp)

    cPickle.dump(np.array(obs_t_verif), open("obs_verif/obs_verif_%s.pkl" % exp_name, 'w'), -1)
    return

if __name__ == "__main__":
    main()
