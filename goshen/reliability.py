
import numpy as np

from util import runConcurrently
from arpsmodelobs import ARPSModelObsFile
from radarobsfile import RadarObsFile

import matplotlib
matplotlib.use('agg')
import pylab

from datetime import datetime, timedelta

def ensembleProbability(ensemble, threshold):
    return (ensemble >= threshold).sum(axis=0).astype(float) / ensemble.shape[0]

def computeReliability(base_path, obs_path, t_ens, base_time, radar, n_ens_members, prob_bins, refl_thresholds):
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
#           print "Loading member %d ..." % (n_ens + 1)
            model_obs = ARPSModelObsFile("%s/%s%03d%s%06d" % (base_path, radar, n_ens + 1, state, t_ens), mpi_config=(2, 12))
            ens_obs.append(model_obs)
    except IOError:
        return np.nan

    ens_refl = np.empty((len(ens_obs),) + ens_obs[0]['Z'].shape, dtype=float)
    for idx in range(len(ens_obs)):
        ens_refl[idx] = ens_obs[idx]['Z']

    radar_obs = RadarObsFile("%s/%s.%s" % (obs_path, radar, (base_time + timedelta(seconds=t_ens)).strftime("%Y%m%d.%H%M%S")))
    obs_refl = radar_obs['Z']

    reliabilities = []
    samples = []

    for r_thresh in refl_thresholds:
        reliability = np.empty((prob_bins.shape[0] - 1,), dtype=float)
        sample = np.empty((prob_bins.shape[0] - 1,), dtype=float)
        ens_prob = ensembleProbability(ens_refl, r_thresh)
        for bin_no in range(prob_bins.shape[0] - 1):
            bin_lb = prob_bins[bin_no]
            bin_ub = prob_bins[bin_no + 1]

            if bin_no < prob_bins.shape[0] - 2:
                bin_points = ((ens_prob >= bin_lb) & (ens_prob < bin_ub))
            else:
                bin_points = ((ens_prob >= bin_lb) & (ens_prob <= bin_ub))

            verif_bin_points = (obs_refl >= r_thresh) & bin_points

            n_points = bin_points.sum()
            n_verif_points = verif_bin_points.sum()

            reliability[bin_no] = float(n_verif_points) / n_points
            sample[bin_no] = n_points

        reliabilities.append(reliability)
        samples.append(sample)

    return reliabilities, samples

def plotReliabilities(reliabilities, samples, bins, thresholds, file_name):
    pylab.figure()

    probs = (bins[1:] + bins[:-1]) / 2
    for thresh, reliability in zip(thresholds, reliabilities):
        pylab.plot(probs, reliability, '-', label=thresh)

    pylab.plot([0., 1.], [0., 1.], 'k--')

    pylab.legend(loc=2)

    pylab.xlim([0, 1])
    pylab.ylim([0, 1])

    pylab.savefig(file_name)
    pylab.close()

def main():
    base_time = datetime(2009, 6, 5, 18, 0, 0)
    exp_names_1km = [ "1kmf-zupdtpt", "1kmf-z-no-05XP", "1kmf-z-no-mm", "1kmf-z-no-snd", "1kmf-z-no-mm-05XP", "1kmf-sndr0h=50km" ]
    exp_names_3km = [ "3km-fixed-radar", "3kmf-r0h=12km", "3kmf-pr0h=16km", "3kmf-r0h=18km", "3km-7dBZ,5ms" ]#, "3kmf-n0r=2e6" ]

    data_path = "/caps2/tsupinie/"
    exp_names = exp_names_1km
    obs_path = "/data6/tsupinie/goshen/qc/1km/"
    t_ens_start = 14400
    t_ens_stop = 18000
    t_ens_step = 300
    n_ens_members = 40

    prob_bins = np.arange(0, 1.1, 0.1)
    refl_thresholds = [ 45 ] #np.arange(20, 50, 5, dtype=int)

    times = range(t_ens_start, t_ens_stop + t_ens_step, t_ens_step)
    reliabilities = []
    samples = []

    for exp_name in exp_names:
        print exp_name
        base_path = "%s%s/" % (data_path, exp_name)
        reliability, sample = runConcurrently(computeReliability, times, max_concurrent=8, zip_result=True, args=(base_path, obs_path, "__placeholder__", base_time, 'KCYS', n_ens_members, prob_bins, refl_thresholds))

        reliabilities.append(reliability)
        samples.append(sample)

    reliabilities = zip(*reliabilities)
    samples = zip(*samples)

    for time, rel, samp in zip(times, reliabilities, samples):
        plotReliabilities([ r[0] for r in rel ], [ s[0] for s in samp ], prob_bins, exp_names, "reliability_%02ddBZ_%06d.png" % (refl_thresholds[0], time))
    return

if __name__ == "__main__":
    main()
