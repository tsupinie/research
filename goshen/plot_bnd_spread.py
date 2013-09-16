
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

import cPickle
from math import ceil

from temporal import goshen_3km_temporal
from util import publicationFigure

def main():
    exp_names = [ "3km-fixed-radar", "3kmf-7dBZ,5ms", "3kmf-r0h=12km", "3kmf-mult=1.03" ]
    bnd_spread = {}

    temp = goshen_3km_temporal(start=10800, end=14400)

    for exp in exp_names:
        bnd_spread[exp] = {}
        bnd_spread[exp]['anal'] = cPickle.load(open("spread_bnd_%s.pkl" % exp, 'r'))

        try:
            bnd_spread[exp]['fcst'] = cPickle.load(open("spread_bnd_%s_fcst.pkl" % exp, 'r'))
        except IOError:
            print "File not found: 'spread_bnd_%s_fcst.pkl'" % exp 
            bnd_spread[exp]['fcst'] = dict( (k, v) for k, v in bnd_spread[exp]['anal'].iteritems() )

    pylab.figure(figsize=(16, 12))
    pylab.subplots_adjust(left=0.07, bottom=0.08, right=0.97, top=0.95, wspace=0.20)

    def subplotFactory(var):
        variables = {'u':r'$u$', 'w':r'$w$', 'pt':r'$\theta$', 'qv':r'$q_v$'}
        def doSubplot(multiplier=1):
            for exp, spd in bnd_spread.iteritems():
                merged_spd = np.vstack((spd['fcst'][var], spd['anal'][var])).flatten('F')
                times = np.array(temp.getTimes())[:, np.newaxis].repeat(2, axis=1).flatten('C')

                pylab.plot(times, merged_spd, label=exp)
            pylab.xlim(temp[0], temp[-1])
            pylab.xticks(temp.getTimes(), temp.getStrings("%H%M", aslist=True), size=(12 * multiplier), rotation=30)
            pylab.xlabel("Time (UTC)", size=(12 * multiplier))

            if var == 'qv':
                pylab.yticks(pylab.yticks()[0],  pylab.yticks()[0] * 1000.)
            pylab.yticks(size=(12 * multiplier))
            pylab.ylabel("Spread in %s" % variables[var], size=(12 * multiplier))
        return doSubplot

    funcs = []
    for var in ['u', 'w', 'pt', 'qv']:
        funcs.append(subplotFactory(var))
    publicationFigure(funcs, (2, 2), corner='ur')

#   pylab.suptitle("Boundary Spread Comparison")
    pylab.savefig("spread_bnd.png")
    return

if __name__ == "__main__":
    main()
