
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

from arpsmodelobs import ARPSModelObsFile
from radarobsfile import RadarObsFile
from util import runConcurrently

from datetime import datetime, timedelta
import os
import cPickle

def ensembleProbability(ensemble, threshold):
    return (ensemble >= threshold).sum(axis=0).astype(float) / ensemble.shape[0]

def POD(ensemble_prob, obs, prob_threshold, refl_threshold, points):
    hits =   ((ensemble_prob >= prob_threshold) & (obs >= refl_threshold)) & points
    misses = ((ensemble_prob < prob_threshold) & (obs >= refl_threshold)) & points

    n_hits = hits.sum()
    n_misses = misses.sum()

#   print "thresh, n_hits, n_misses =", threshold, n_hits, n_misses
#   print "POD =", float(n_hits) / (n_hits + n_misses)

    return float(n_hits) / (n_hits + n_misses)

def POFD(ensemble_prob, obs, prob_threshold, refl_threshold, points):
    false_alarms = ((ensemble_prob >= prob_threshold) & (obs < refl_threshold)) & points
    correct_negs = ((ensemble_prob < prob_threshold) & (obs < refl_threshold)) & points

    n_false_alarms = false_alarms.sum()
    n_correct_negs = correct_negs.sum()

#   print "thresh, n_false_alarms, n_correct_negs =", threshold, n_false_alarms, n_correct_negs
#   print "POFD =", float(n_false_alarms) / (n_correct_negs + n_false_alarms)

    return float(n_false_alarms) / (n_correct_negs + n_false_alarms)

def findROCs(base_path, obs_path, t_ens, base_time, radar, n_ens_members, prob_thresholds, refl_thresholds):
    ens_obs = []
    state = "an"
    print "Working on time %06d ..." % t_ens
    try:
        for n_ens in range(n_ens_members):
#           print "Loading member %d ..." % (n_ens + 1)
            model_obs = ARPSModelObsFile("%s/%s%03d%s%06d" % (base_path, radar, n_ens + 1, state, t_ens))
            ens_obs.append(model_obs)
    except AssertionError:
        for n_ens in range(n_ens_members):
            print "Loading member %d ..." % (n_ens + 1)
            model_obs = ARPSModelObsFile("%s/%s%03d%s%06d" % (base_path, radar, n_ens + 1, state, t_ens), mpi_config=(2, 12))
            ens_obs.append(model_obs)
    except IOError as e:
        print "IOError:", e
        return [ [ (np.nan, np.nan), (np.nan, np.nan) ] ]

    ens_refl = np.empty((len(ens_obs),) + ens_obs[0]['Z'].shape, dtype=float)
    for idx in range(len(ens_obs)):
        ens_refl[idx] = ens_obs[idx]['Z']

    rof_name = "%s/%s.%s" % (obs_path, radar, (base_time + timedelta(seconds=t_ens)).strftime("%Y%m%d.%H%M%S"))
    radar_obs = RadarObsFile(rof_name)
    obs_refl = radar_obs['Z']

    good_obs = (obs_refl >= 15) & np.any(ens_refl >= 15, axis=0) #& ((obs_refl[np.newaxis, ...] >= 15) | (ens_refl >= 15))

    rocs = []

    for r_thresh in refl_thresholds:
        ens_prob = ensembleProbability(ens_refl, r_thresh)

        roc = [(1., 1.)]
    
        for p_thresh in prob_thresholds:
            pod = POD(ens_prob, obs_refl, p_thresh, r_thresh, good_obs)
            pofd = POFD(ens_prob, obs_refl, p_thresh, r_thresh, good_obs)

            if np.isnan(pod): pod = 0.
            if np.isnan(pofd): pofd = 1.

            roc.append((pod, pofd))
        roc.append((0., 0.))
        rocs.append(roc)
    return rocs

def plotROCs(rocs, prob_thresholds, refl_thresholds, file_name):
    if type(rocs[0]) == tuple:
        rocs = [ rocs ]

    matplotlib.rcParams['axes.color_cycle'] = ['k', 'r', 'g', 'b', 'm', 'c']
    pylab.figure()

    for roc, thresh in zip(rocs, refl_thresholds):
        pod, pofd = zip(*roc)
        pylab.plot(pofd, pod, label=thresh)

#   pylab.plot(pofd, pod, 'ko')
#   for thresh, (pod, pofd) in zip(thresholds, roc):
#       pylab.text(pofd + 0.02, pod - 0.02, "%d dBZ" % thresh, ha='left', va='top')

    pylab.plot([0., 1.], [0., 1.], 'k--')

    pylab.xlim([0, 1])
    pylab.ylim([0, 1])

    pylab.legend(loc=4)

    pylab.savefig(file_name)
    pylab.close()
    return

def computeAUC(roc):
    pod, pofd = zip(*roc)
    return -np.trapz(pod, x=pofd)

def plotAUCs(AUCs, times, exp_names, base_time, title, file_name, exclude=[]):
    pylab.figure()
    pylab.axes((0.1, 0.125, 0.8, 0.8))

    for AUC, exp in zip(AUCs, exp_names):
        if exp not in exclude:
            pylab.plot(times, AUC, label=exp)

    pylab.axhline(y=0.5, color='k', linestyle=':', linewidth=0.5)

    pylab.legend(loc=3)

    pylab.xlim([times[0], times[-1]])
    pylab.ylim([0, 1])

    pylab.xticks(times, [ (base_time + timedelta(seconds=int(t))).strftime("%H%M") for t in times ], rotation=30, size='large')
    pylab.yticks(size='large')

    pylab.xlabel("Time (UTC)", size='large')
    pylab.ylabel("AUC", size='large')

    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def main():
    base_time = datetime(2009, 6, 5, 18, 0, 0)
    exp_names_1km = [ "1kmf-sndr0h=25km", "1kmf-zs25-no-05XP", "1kmf-z-no-snd", "1kmf-zs25-no-mm", "1kmf-zs25-no-mm-05XP", "1kmf-z-no-v2" ] # "1kmf-newbc", "1kmf-r0h=4km", "1kmf-bc7dBZ,5ms", "1kmf-bcmult=1.03", "1kmf-r0h=6km-bc7dBZ,5ms", "1kmf-5dBZ,3ms-bc7dBZ,5ms"
    exp_names_3km = [ "3km-fixed-radar", "3kmf-r0h=12km", "3kmf-pr0h=16km", "3kmf-r0h=18km", "3km-7dBZ,5ms" ]#, "3kmf-n0r=2e6" ]
    plot_exp_names_1km = { '1kmf-control':"Control BC", '1kmf-newbc':"Control", '1kmf-r0h=4km':"r$_{0h}$=4km", '1kmf-prtrgn=1':"Storm Pert.", '1kmf-sndr0h=25km':"CTRL", '1kmf-zs25-no-05XP':"NO_MWR", '1kmf-zs25-no-mm-05XP':"NO_MWR_MM", '1kmf-zs25-no-mm':"NO_MM", '1kmf-z-no-snd':"NO_SND", '1kmf-z-no-v2':"NO_V2 ", '1kmf-sndr0h=50km':"Sounding $r_{0h}$ = 50 km", '1kmf-bc7dBZ,5ms':"r$_{0h}$=4km, BC: 7dBZ,5ms", '1kmf-bcmult=1.03':"BC: mult=1.03", '1kmf-r0h=6km-bc7dBZ,5ms':"BC: 7dBZ,5ms", '1kmf-5dBZ,3ms-bc7dBZ,5ms':"5dBZ,3ms, BC: 7dBZ,5ms" }
    plot_exp_names_3km = { '3km-fixed-radar':"Control", '3kmf-r0h=12km':"r$_{0h}$=12km", '3kmf-r0h=18km':"r$_{0h}$=8km", '3kmf-pr0h=16km':"Init r$_{0h}$=16km", '3km-7dBZ,5ms':"7 dBZ, 5 m s$^{-1}$"}

    data_path = "/caps2/tsupinie/"
    exp_names = exp_names_1km
    plot_exp_names = plot_exp_names_1km
    obs_path = "/data6/tsupinie/goshen/qc/1km/"
    t_ens_start = 14400
    t_ens_stop = 18000
    t_ens_step = 300
    n_ens_members = 40

    overwrite = False

    prob_thresholds = np.arange(0, 1.01, 0.01)
    refl_thresholds = [ 45 ] # np.arange(20, 50, 5, dtype=int)

    times = range(t_ens_start, t_ens_stop + t_ens_step, t_ens_step)

    rocs = []
    for exp_name in exp_names:
        print exp_name
        pkl_name = "roc_pkl/%s_%02ddBZ_roc.pkl" % (exp_name, refl_thresholds[0])
        if not overwrite and os.path.exists(pkl_name):
            roc = cPickle.load(open(pkl_name, 'r'))
        else:
            base_path = "%s%s/" % (data_path, exp_name)
            roc = runConcurrently(findROCs, times, max_concurrent=8, args=(base_path, obs_path, '__placeholder__', base_time, 'KCYS', n_ens_members, prob_thresholds, refl_thresholds))
            cPickle.dump(roc, open(pkl_name, 'w'), -1)
        rocs.append(roc)

    rocs = zip(*rocs)
    print len(rocs), len(times)
    AUCs = []
    for time, roc_group in zip(times, rocs):
        AUC = [ computeAUC(r[0]) for r in roc_group ]
        plotROCs([ r[0] for r in roc_group ], prob_thresholds, exp_names, "roc_1km_%02ddBZ_%06d.png" % (refl_thresholds[0], time))

        print time, AUC
        AUCs.append(AUC)

#   plotAUCs(zip(*AUCs), times, [ plot_exp_names[e] for e in exp_names], base_time, "ROC Area for $Z$ > %02d dBZ" % refl_thresholds[0], "roc_1km_%02ddBZ_AUC.png" % refl_thresholds[0])
#   plotAUCs(zip(*AUCs), times, [ plot_exp_names[e] for e in exp_names], base_time, "ROC Area for $Z$ > %02d dBZ" % refl_thresholds[0], "roc_1km_%02ddBZ_AUC_opt.png" % refl_thresholds[0], 
#       exclude=[ plot_exp_names[e] for e in ['1kmf-no-05XP', '1kmf-no-mm-05XP'] ])

#   plot_exp_names['1kmf-r0h=4km'] = "All V2 Data"

    plotAUCs(zip(*AUCs), times, [ plot_exp_names[e] for e in exp_names], base_time, "ROC Area for $Z$ > %02d dBZ" % refl_thresholds[0], "roc_1km_%02ddBZ_AUC_dd.png" % refl_thresholds[0], 
        exclude=[ ])
    return

if __name__ == "__main__":
    main()
