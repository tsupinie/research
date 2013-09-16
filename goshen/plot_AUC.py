
import numpy as np

from roc import computeAUC
from temporal import goshen_1km_temporal
from util import publicationFigure

import cPickle
from collections import OrderedDict

import matplotlib
#matplotlib.use('agg')
import pylab

#   plotAUCs(zip(*AUCs), times, [ plot_exp_names[e] for e in exp_names], base_time, "ROC Area for $Z$ > %02d dBZ" % refl_thresholds[0], "roc_1km_%02ddBZ_AUC_dd.png" % refl_thresholds[0], 
def main():
    temp = goshen_1km_temporal(start=14400)
    experiments = OrderedDict([('1kmf-sndr0h=25km', 'CTRL'), ('1kmf-zs25-no-05XP', 'NO_MWR'), ('1kmf-z-no-snd', 'NO_SND'), ('1kmf-zs25-no-mm', 'NO_MM'), ('1kmf-zs25-no-mm-05XP', 'NO_MWR_MM'), ('1kmf-z-no-v2', 'NO_V2')])
    refl_thresh = [ 25, 45 ]

    all_AUCs = []
    for thresh in refl_thresh:
        rocs = []
        for exp in experiments.iterkeys():
            pkl_name = "roc_pkl/%s_%02ddBZ_roc.pkl" % (exp, thresh)
            roc = cPickle.load(open(pkl_name, 'r'))
            rocs.append(roc)

        rocs = zip(*rocs)
        AUCs = []
        for time, roc_group in zip(temp, rocs):
            AUC = [ computeAUC(r[0]) for r in roc_group ]
            AUCs.append(AUC)
        all_AUCs.append(AUCs)

    all_AUCs = [ zip(*a) for a in all_AUCs ]

    def subplotFactory(AUC_group, thresh):
        colors = dict(zip(experiments.iterkeys(), ['k', 'r', 'g', 'b', 'c', 'm']))
        def doSubplot(multiplier=1.0, layout=(-1, -1)):
            for idx, (exp, exp_name) in enumerate(experiments.iteritems()):
                pylab.plot(temp.getTimes(), AUC_group[idx], label=exp_name, color=colors[exp])

            pylab.axhline(y=0.5, color='k', linestyle=':')

            n_row, n_col = layout
            if n_row == 2:
                pylab.xlabel("Time (UTC)", size='large')
                pylab.xlim(temp.getTimes()[0], temp.getTimes()[-1])
                pylab.xticks(temp.getTimes(), temp.getStrings("%H%M", aslist=True), rotation=30, size='x-large')

                pylab.legend(loc=3)

            else:
                pylab.xlim(temp.getTimes()[0], temp.getTimes()[-1])
                pylab.xticks(temp.getTimes(), [ "" for t in temp ])

            pylab.ylabel(r"AUC (Reflectivity $\geq$ %d dBZ)" % thresh, size='large')
            pylab.ylim(0.4, 1.0)
            pylab.yticks(size='x-large')
            return

        return doSubplot

    subplots = []
    for AUC_group, thresh in zip(all_AUCs, refl_thresh):
        subplots.append(subplotFactory(AUC_group, thresh))

    pylab.figure(figsize=(8, 9.5))
    pylab.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1, hspace=0.05)
    publicationFigure(subplots, (2, 1), corner='ur')

    pylab.savefig("roc_AUC.png")
    pylab.close()
    return

if __name__ == "__main__":
    main()
