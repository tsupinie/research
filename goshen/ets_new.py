
from multiprocessing import Process, Queue
from datetime import datetime, timedelta
import cPickle

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib.ticker import FixedFormatter, FixedLocator
from matplotlib.colors import LinearSegmentedColormap

from arpsmodelobs import ARPSModelObsFile
from radarobsfile import RadarObsFile
from grid import goshen_3km_grid

def plotConfusion(confusion, grid, title, file_name, inset=None, fudge=16):
    pylab.figure()
    axmain = pylab.axes()

    confusion_cmap_dict = {
        'red':(
            (0.0, 1.0, 1.0), 
            (0.2, 1.0, 0.5), # Missing to Correct Negative
            (0.4, 0.5, 1.0), # Correct Negative to False Alarm
            (0.6, 1.0, 1.0), # False Alarm to Miss
            (0.8, 1.0, 0.0), # Miss to Hit
            (1.0, 0.0, 0.0),
        ),
        'green':(
            (0.0, 1.0, 1.0),
            (0.2, 1.0, 0.5), 
            (0.4, 0.5, 0.0),
            (0.6, 0.0, 1.0),
            (0.8, 1.0, 0.5),
            (1.0, 0.5, 0.5),
        ),
        'blue':(
            (0.0, 1.0, 1.0),
            (0.2, 1.0, 0.5), 
            (0.4, 0.5, 1.0),
            (0.6, 1.0, 0.0),
            (0.8, 0.0, 0.0),
            (1.0, 0.0, 0.0),
        ),
    }
    confusion_cmap = LinearSegmentedColormap('confusion', confusion_cmap_dict, 256)

    xs, ys = grid.getXY()

    pylab.pcolormesh(xs, ys, confusion, cmap=confusion_cmap, vmin=-1, vmax=3)

    tick_labels = [ "Missing", "Correct\nNegative", "False\nAlarm", "Miss", "Hit" ]
    tick_locs = np.linspace(-1, len(tick_labels) - 3, len(tick_labels))
    tick_locs += (tick_locs[1] - tick_locs[0]) / 2
    bar = pylab.colorbar()
    bar.locator = FixedLocator(tick_locs)
    bar.formatter = FixedFormatter(tick_labels)
    pylab.setp(pylab.getp(bar.ax, 'ymajorticklabels'), fontsize='large')
    bar.update_ticks()

    grid.drawPolitical()

    if inset:
        lb_y, lb_x = [ b.start for b in inset ]
        ub_y, ub_x = [ b.stop + fudge for b in inset ]

        inset_exp = (slice(lb_y, ub_y), slice(lb_x, ub_x))

        axins = zoomed_inset_axes(pylab.gca(), 2, loc=4)
        pylab.sca(axins)

        pylab.pcolormesh(xs[inset_exp], ys[inset_exp], confusion[inset_exp], cmap=confusion_cmap, vmin=0, vmax=3)
        grid.drawPolitical()

        pylab.xlim([lb_x * gs_x, (ub_x - 1) * gs_x])
        pylab.ylim([lb_y * gs_y, (ub_y - 1) * gs_y])

        mark_inset(axmain, axins, loc1=1, loc2=3, fc='none', ec='k')

    pylab.sca(axmain)
    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return
    
def ETS(forecast, observation, threshold, good_markers):
    hits = (forecast >= threshold) & (observation >= threshold) & good_markers
    misses = (forecast < threshold) & (observation >= threshold) & good_markers
    false_alarms = (forecast >= threshold) & (observation < threshold) & good_markers
    correct_negatives = (forecast < threshold) & (observation < threshold) & good_markers

    num_hits = hits.sum()
    num_misses = misses.sum()
    num_fa = false_alarms.sum()
    num_cn = correct_negatives.sum()

    total = num_hits + num_misses + num_fa + num_cn

    hits_random = float(num_hits + num_misses) * (num_hits + num_fa) / total
    ets = float(num_hits - hits_random) / (num_hits + num_misses + num_fa - hits_random)

    confusion = -np.ones(observation.shape, dtype=np.int32)

    confusion = np.where(correct_negatives, 0, confusion)
    confusion = np.where(false_alarms, 1, confusion)
    confusion = np.where(misses, 2, confusion)
    confusion = np.where(hits, 3, confusion)

    return ets, confusion

def doETS(radar, model_path, obs_path, t_ens, base_time, refl_threshold, vel_threshold, grid, pipe):
    try:
        model_obs = ARPSModelObsFile("%s/%san%06d" % (model_path, radar, t_ens))
    except AssertionError:
        model_obs = ARPSModelObsFile("%s/%san%06d" % (model_path, radar, t_ens), mpi_config=(2, 12))
    except IOError:
        pipe.put((t_ens, np.nan, np.nan), block=False)
        return

    radar_obs = RadarObsFile("%s/%s.%s" % (obs_path, radar, (base_time + timedelta(seconds=t_ens)).strftime("%Y%m%d.%H%M%S")))
    obs_refl = np.transpose(radar_obs.reflectivity, [2, 0, 1])
    obs_vel = np.transpose(radar_obs.radial_velocity, [2, 0, 1])
    good_markers_refl = (model_obs['vr'] > -90)
    good_markers_vel = (model_obs['Z'] > 15) & (obs_refl > 15)# & ((model_obs['Z'] > 15) | (obs_refl > 15))

    ets_refl, confusion_refl = ETS(model_obs['Z'], obs_refl, refl_threshold, good_markers_refl)
    ets_vel, confusion_vel = ETS(np.abs(model_obs['vr']), np.abs(obs_vel), vel_threshold, good_markers_vel)

#   print np.isnan(confusion_refl).sum()

#   if t_ens == 14400 and radar == "KCYS":
#       xs, ys = grid.getXY()
#       for tilt in range(model_obs['Z'].shape[0]):
#           plotConfusion(confusion_refl[tilt], grid, r"Confusion for $Z$ > %02d dBZ" % refl_threshold, "confusion_tilt%02d_%02ddBZ.png" % (tilt, refl_threshold))
#           plotConfusion(confusion_vel[tilt], grid, r"Confusion for |$v_r$| > %02d m s$^{-1}$" % vel_threshold, "confusion_tilt%02d_%02dms.png" % (tilt, vel_threshold))
#
#           pylab.figure()
#           pylab.pcolormesh(xs, ys, np.where(good_markers_vel[tilt], np.abs(model_obs['vr'][tilt]), np.nan), cmap=matplotlib.cm.get_cmap('winter'), vmin=0, vmax=40)
#           pylab.colorbar()
#           pylab.savefig("vr_model_tilt_%02d.png" % tilt)
#           pylab.close()
#
#           pylab.figure()
#           pylab.pcolormesh(xs, ys, np.where(good_markers_vel[tilt], np.abs(obs_vel[tilt]), np.nan), cmap=matplotlib.cm.get_cmap('winter'), vmin=0, vmax=40)
#           pylab.colorbar()
#           pylab.savefig("vr_obs_tilt_%02d.png" % tilt)
#           pylab.close()

    pipe.put((t_ens, ets_refl, ets_vel), block=False)
    return

def main():
    model_paths = [ "/caps1/tsupinie/3km-fixed-radar/", "/caps2/tsupinie/3km-control/", "/caps2/tsupinie/3km-n0r=8e5/", "/caps2/tsupinie/3km-7dBZ,5ms/" ]
    obs_path = "/data6/tsupinie/goshen/qc/3km/"
    t_ens_start = 14400
    t_ens_stop = 18000
    t_ens_step = 300

    refl_threshold = 40
    vel_threshold = 30

    base_time = datetime(2009, 6, 5, 18, 0, 0)

    grid = goshen_3km_grid()

    all_ets_refl = {}
    all_ets_vel = {}
    for model_path in model_paths:
        exp_key = model_path.split("/")[-2]
        all_ets_refl[exp_key] = {}
        all_ets_vel[exp_key] = {}
        for radar in [ 'KCYS', 'KFTG', 'KRIW' ]:
            procs = []
            ets = []
            pipe = Queue(len(range(t_ens_start, t_ens_stop + t_ens_step, t_ens_step)))
            for t_ens in xrange(t_ens_start, t_ens_stop + t_ens_step, t_ens_step):
                proc = Process(target=doETS, name="%s at %06d" % (radar, t_ens), args=(radar, model_path, obs_path, t_ens, base_time, refl_threshold, vel_threshold, grid, pipe))
                proc.start()
                procs.append(proc)

            proc_names = [ p.name for p in procs ]         

            while len(procs) > 0:
                ets_time = pipe.get()
                time, ets_refl, ets_vel = ets_time

                proc_idx = proc_names.index("%s at %06d" % (radar, time))
                procs.pop(proc_idx)
                proc_names.pop(proc_idx)

                ets.append(ets_time)

            all_ets_refl[exp_key][radar] = zip(*sorted(ets, key=lambda x: x[0]))[1]
            all_ets_vel[exp_key][radar] = zip(*sorted(ets, key=lambda x: x[0]))[2]

    cPickle.dump(all_ets_refl, open("all_ets_new_%02ddBZ.pkl" % refl_threshold, 'w'), -1)
    cPickle.dump(all_ets_vel, open("all_ets_new_%02dms.pkl" % vel_threshold, 'w'), -1)
    return

if __name__ == "__main__":
    main()
