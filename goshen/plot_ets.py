
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

import cPickle
from datetime import datetime, timedelta

def main():
    exp_names = ['mod-05XP', 'mm', 'no-mm']
    labels = { 'mod-05XP':"MM + MWR05XP", 'mm':'MM', 'no-mm':'No MM' }
    threshold = 40

    times = np.arange(14400, 18300, 300)
    base_time = datetime(2009, 6, 5, 18, 0, 0)

    time_start = np.where(times == 14400)[0][0]

    pylab.figure(figsize=(8, 6))
    pylab.axes((0.1, 0.125, 0.8, 0.8))
    for name in exp_names:
        mean_ets = cPickle.load(open("%s_%ddBZ.pkl" % (name, threshold), 'r'))

        pylab.plot(times[time_start:], mean_ets[time_start:], label=labels[name])

    y_lb, y_ub = pylab.ylim(0, 0.7)
    pylab.plot([14400, 14400], [y_lb, y_ub], 'k--')
    pylab.ylim(y_lb, y_ub)
    pylab.xlim(times[time_start], times[-1])

    pylab.legend(loc=2)
    pylab.suptitle("ETS for Probability Matched Mean of Reflectivity (%d dBZ threshold)" % threshold)

    pylab.xlabel("Time (UTC)", size='large')
    pylab.ylabel("ETS", size='large')
    pylab.xticks(times[time_start:], [ (base_time + timedelta(seconds=int(t))).strftime("%H%M") for t in times[time_start:] ], rotation=30, size='large')
    pylab.yticks(size='large')

    pylab.savefig("all_pmmean_ets_%ddBZ.png" % threshold)
    return

if __name__ == "__main__":
    main()
