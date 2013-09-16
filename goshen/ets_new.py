
from datetime import datetime, timedelta
import cPickle
import glob

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib.ticker import FixedFormatter, FixedLocator
from matplotlib.colors import LinearSegmentedColormap

from util import runConcurrently
from arpsmodelobs import ARPSModelObsFile
from radarobsfile import RadarObsFile
from grid import goshen_1km_grid, goshen_3km_grid
from temporal import goshen_1km_temporal

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
    tick_labels = [ "Missing", "Correct\nNegative", "False\nAlarm", "Miss", "Hit" ]
    min_label = -1

    xs, ys = grid.getXY()

    pylab.pcolormesh(xs, ys, confusion, cmap=confusion_cmap, vmin=min_label, vmax=(min_label + len(tick_labels) - 1))

    tick_locs = np.linspace(-1, min_label + len(tick_labels) - 2, len(tick_labels))
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
    
def ETS(forecast, observation, threshold, good_markers, bbox=(slice(None), slice(None), slice(None))):
    hits = (forecast >= threshold) & (observation >= threshold) & good_markers
    misses = (forecast < threshold) & (observation >= threshold) & good_markers
    false_alarms = (forecast >= threshold) & (observation < threshold) & good_markers
    correct_negatives = (forecast < threshold) & (observation < threshold) & good_markers

    num_hits = hits[bbox].sum()
    num_misses = misses[bbox].sum()
    num_fa = false_alarms[bbox].sum()
    num_cn = correct_negatives[bbox].sum()

    total = num_hits + num_misses + num_fa + num_cn

    hits_random = float(num_hits + num_misses) * (num_hits + num_fa) / total
    ets = float(num_hits - hits_random) / (num_hits + num_misses + num_fa - hits_random)

    confusion = -np.ones(observation.shape, dtype=np.int32)

    confusion = np.where(correct_negatives, 0, confusion)
    confusion = np.where(false_alarms, 1, confusion)
    confusion = np.where(misses, 2, confusion)
    confusion = np.where(hits, 3, confusion)

    return ets, confusion

def ETSens(forecasts, observation, prob_threshold, refl_threshold, good_markers, bbox=(slice(None), slice(None), slice(None))):
    ens_prob = (forecasts >= refl_threshold).sum(axis=0) / float(len(forecasts))

    hits = (ens_prob >= prob_threshold) & (observation >= refl_threshold) & good_markers
    misses = (ens_prob < prob_threshold) & (observation >= refl_threshold) & good_markers
    false_alarms = (ens_prob >= prob_threshold) & (observation < refl_threshold) & good_markers
    correct_negatives = (ens_prob < prob_threshold) & (observation < refl_threshold) & good_markers

    num_hits = hits[bbox].sum()
    num_misses = misses[bbox].sum()
    num_fa = false_alarms[bbox].sum()
    num_cn = correct_negatives[bbox].sum()

    total = num_hits + num_misses + num_fa + num_cn

    hits_random = float(num_hits + num_misses) * (num_hits + num_fa) / total
    ets = float(num_hits - hits_random) / (num_hits + num_misses + num_fa - hits_random)

    confusion = -np.ones(observation.shape, dtype=np.int32)

    confusion = np.where(correct_negatives, 0, confusion)
    confusion = np.where(false_alarms, 1, confusion)
    confusion = np.where(misses, 2, confusion)
    confusion = np.where(hits, 3, confusion)

    return ets, confusion

def doETS(radar, model_path, obs_path, t_ens, base_time, refl_threshold, vel_threshold, grid, n_ens_members=-1, prob_threshold=0.5, bbox=(slice(None), slice(None), slice(None))):
    try:
        radar_obs = RadarObsFile("%s/%s.%s" % (obs_path, radar, (base_time + timedelta(seconds=t_ens)).strftime("%Y%m%d.%H%M%S")))
    except IOError:
        return np.nan, np.nan, np.nan, np.nan

    if n_ens_members <= 0:
        try:
            model_obs = ARPSModelObsFile("%s/%san%06d" % (model_path, radar, t_ens))
        except AssertionError:
            model_obs = ARPSModelObsFile("%s/%san%06d" % (model_path, radar, t_ens), mpi_config=(2, 12))
        except IOError:
            return np.nan, np.nan, np.nan, np.nan

        good_markers_refl = (model_obs['vr'] > -90)
        good_markers_vel = (model_obs['Z'] > 15) & (radar_obs['Z'] > 15)# & ((model_obs['Z'] > 15) | (obs_refl > 15))

        ets_refl, confusion_refl = ETS(model_obs['Z'], radar_obs['Z'], refl_threshold, good_markers_refl, bbox=bbox)
        ets_vel, confusion_vel = ETS(np.abs(model_obs['vr']), np.abs(radar_obs['vr']), vel_threshold, good_markers_vel, bbox=bbox)
    else:
        model_obs = []
        for n_ens in range(n_ens_members):
            try:
                model_obs.append(ARPSModelObsFile("%s/%s%03dan%06d" % (model_path, radar, n_ens + 1, t_ens)))
            except AssertionError:
                model_obs.append(ARPSModelObsFile("%s/%s%03dan%06d" % (model_path, radar, n_ens + 1, t_ens), mpi_config=(2, 12)))
            except IOError:
                return np.nan, np.nan, np.nan, np.nan

        model_refl = np.array([ m['Z'] for m in model_obs ])
        model_vel = np.array([ m['vr'] for m in model_obs ])

        good_markers_refl = (model_vel.mean(axis=0) > -90)
        good_markers_vel = (model_refl.mean(axis=0) > 15) & (radar_obs['Z'] > 15)# & ((model_obs['Z'] > 15) | (obs_refl > 15))

        ets_refl, confusion_refl = ETSens(model_refl, radar_obs['Z'], prob_threshold, refl_threshold, good_markers_refl, bbox=bbox)
        ets_vel, confusion_vel = ETSens(np.abs(model_vel), np.abs(radar_obs['vr']), prob_threshold, vel_threshold, good_markers_vel, bbox=bbox)

#       if radar == "KCYS":
#           plotConfusion(confusion_refl[0], grid, "Contingency for %s" % radar, "confusion_%s_tilt00_%06d.png" % (radar, t_ens))
        
    return ets_refl, ets_vel, confusion_refl, confusion_vel

def main():
    base_path = "/caps2/tsupinie/"
    model_paths_3km = [ "/caps1/tsupinie/3km-fixed-radar/", "/caps2/tsupinie/3km-control/", "/caps2/tsupinie/3km-n0r=8e5/", 
        "/caps2/tsupinie/3km-7dBZ,5ms/", "/caps1/tsupinie/3kmf-r0h=12km/", "/caps2/tsupinie/3kmf-pr0h=16km/", "/caps2/tsupinie/3kmf-r0h=18km/" ]
    model_paths_1km = [ "/caps2/tsupinie/1kmf-z-no-snd/"]
    exp_names_1km = [ "1kmf-zupdtpt", "1kmf-z-no-05XP", "1kmf-z-no-mm-05XP", "1kmf-z-no-mm", "1kmf-z-no-v2", "1kmf-z-no-snd" ]

    model_paths = [ "%s%s/" % (base_path, e) for e in exp_names_1km ]
    obs_path = "/data6/tsupinie/goshen/qc/1km/"
    temp = goshen_1km_temporal(start=14400)

    refl_threshold = 40
    vel_threshold = 30

    base_time = datetime(2009, 6, 5, 18, 0, 0)

    grid = goshen_1km_grid()

    all_ets_refl = {}
    all_ets_vel = {}

    confusion_refl = {}
    confusion_vel = {}

    bbox_files = glob.glob("bbox*.pkl")
    bboxes = {}
    bbox_buffer = 10
    bbox_offsets = [0, 10, 0]
    for bfile in bbox_files:
        root, ext = bfile.split(".")
        bits = root.split("_")
        key = "1kmf-%s" % bits[-1]

        bboxes[key] = cPickle.load(open(bfile, 'r'))
        bboxes[key] = (slice(0, 14),) + tuple( slice(b.start - bbox_buffer + o, b.stop + bbox_buffer + o) for b, o in zip(bboxes[key][1:], bbox_offsets[1:]) )

    for model_path in model_paths:
        print model_path
        exp_key = model_path.split("/")[-2]
        all_ets_refl[exp_key] = {}
        all_ets_vel[exp_key] = {}

        confusion_refl[exp_key] = {}
        confusion_vel[exp_key] = {}

        for radar in [ 'KCYS', 'KFTG', '05XP' ]:
            print "%s ..." % radar
            ets_refl, ets_vel, conf_refl, conf_vel = runConcurrently(doETS, temp.getTimes(), zip_result=True, args=(radar, model_path, obs_path, "__placeholder__", base_time, refl_threshold, vel_threshold, grid), 
                kwargs={'n_ens_members':40, 'prob_threshold':0.5, 'bbox':bboxes[exp_key]})

            all_ets_refl[exp_key][radar] = ets_refl
            all_ets_vel[exp_key][radar] = ets_vel

            confusion_refl[exp_key][radar] = conf_refl
            confusion_vel[exp_key][radar] = conf_vel

    cPickle.dump(all_ets_refl, open("all_ets_new_1km_%02ddBZ_bnd.pkl" % refl_threshold, 'w'), -1)
    cPickle.dump(all_ets_vel, open("all_ets_new_1km_%02dms_bnd.pkl" % vel_threshold, 'w'), -1)
    cPickle.dump(confusion_refl, open("all_confusion_1km_%02ddBZ.pkl" % refl_threshold, 'w'), -1)
    cPickle.dump(confusion_vel, open("all_confusion_1km_%02dms.pkl" % vel_threshold, 'w'), -1)
    return

if __name__ == "__main__":
    main()
