
import matplotlib
matplotlib.use('agg')
import pylab

import numpy as np

from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from util import publicationFigure

import cPickle
from collections import OrderedDict

def main():
    experiments = OrderedDict([('1kmf-sndr0h=25km', 'CTRL'), ('1kmf-zs25-no-05XP', 'NO_MWR'), ('1kmf-z-no-snd', 'NO_SND'), ('1kmf-zs25-no-mm', 'NO_MM')]) #, ('1kmf-zs25-no-mm-05XP', 'NO_MWR_MM'), ('1kmf-z-no-v2', 'NO_V2')])

    domain_bounds = (slice(115, 165), slice(105, 155))
    grid = goshen_1km_grid(bounds=domain_bounds)
    domain_bounds = grid.getBounds()

    temporal = goshen_1km_temporal(start=14400)

    tornado_track = zip(*((41.63, -104.383), (41.6134, -104.224)))
    track_xs, track_ys = grid(*[ np.array(list(t)) for t in reversed(tornado_track) ])

    def subplotFactory(exp):
        def doSubplot(multiplier=1.0, layout=(-1, -1)):
            data = cPickle.load(open("updraft_hel/NEPcirc_UH0-3_75m2s2_%s.pkl" % exp, 'r'))

            cmap = matplotlib.cm.get_cmap('RdYlBu_r')
            cmap.set_under('#ffffff')

            xs, ys = grid.getXY()
            pylab.pcolormesh(xs, ys, data[domain_bounds], cmap=cmap, vmin=0.1, vmax=1.0)
            grid.drawPolitical(scale_len=10)

            pylab.plot(track_xs, track_ys, 'mv-', lw=2.5, mfc='k', ms=8)

            pylab.text(0.05, 0.95, experiments[exp], ha='left', va='top', transform=pylab.gca().transAxes, size=18 * multiplier)
            return
        return doSubplot

    subplots = []
    for exp in experiments.iterkeys():
        subplots.append(subplotFactory(exp))

#   pylab.figure(figsize=(12, 8))
    pylab.figure(figsize=(8.5, 8))
    pylab.subplots_adjust(left=0.025, bottom=0.1, right=0.875, top=0.975, hspace=0.1, wspace=0.1)

    publicationFigure(subplots, (2, 2), corner='ur')

    cax = pylab.axes((0.90, 0.1125, 0.020, 0.825))
    bar = pylab.colorbar(cax=cax)
    bar.ax.text(3.5, 0.5, r"Probability of UH $\geq$ 75 m$^2$ s$^{-2}$", rotation=90, transform=bar.ax.transAxes, size='large', va='center')
    bar.ax.tick_params(labelsize='large')

    pylab.savefig("UH0-3_75m2s2_four.png")
    return

if __name__ == "__main__":
    main()
