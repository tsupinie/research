
import matplotlib
matplotlib.use('agg')
import pylab

import numpy as np

from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from util import publicationFigure
from arpsmodelobs import ARPSModelObsFile

import cPickle
import argparse
from collections import OrderedDict

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ens-time', dest='ens_time', required=True, type=int)

    args = ap.parse_args()

    base_path="/caps2/tsupinie/"
    experiments = OrderedDict([('1kmf-sndr0h=25km', 'CTRL'), ('1kmf-zs25-no-05XP', 'NO_MWR'), ('1kmf-z-no-snd', 'NO_SND'), ('1kmf-zs25-no-mm', 'NO_MM'), ('1kmf-zs25-no-mm-05XP', 'NO_MWR_MM'), ('1kmf-z-no-v2', 'NO_V2')])
    domain_bounds = (slice(105, 155), slice(110, 160))
    grid = goshen_1km_grid(bounds=domain_bounds)
    domain_bounds = grid.getBounds()
    temporal = goshen_1km_temporal(start=14400)

    tornado_track = zip(*((41.63, -104.383), (41.6134, -104.224)))
    track_xs, track_ys = grid(*[ np.array(list(t)) for t in reversed(tornado_track) ])
    t_ens = args.ens_time

    def subplotFactory(exp, time_sec):
        def doSubplot(multiplier=1.0, layout=(-1, -1)):
            data = cPickle.load(open("cold_pool_%s.pkl" % exp, 'r'))
            wdt = temporal.getTimes().index(time_sec)

            try:
                mo = ARPSModelObsFile("%s/%s/KCYSan%06d" % (base_path, exp, time_sec))
            except AssertionError:
                mo = ARPSModelObsFile("%s/%s/KCYSan%06d" % (base_path, exp, time_sec), mpi_config=(2, 12))
            except:
                print "Can't load reflectivity ..."
                mo = {'Z':np.zeros((1, 255, 255), dtype=np.float32)}

            cmap = matplotlib.cm.get_cmap('Blues_r')
            cmap.set_over('#ffffff')
#           cmap.set_under(tuple( cmap._segmentdata[c][0][-1] for c in ['red', 'green', 'blue'] ))
#           cmap.set_under(cmap[0])

            xs, ys = grid.getXY()
#           pylab.pcolormesh(xs, ys, data['t'][2][domain_bounds], cmap=cmap, vmin=288., vmax=295.)
            pylab.contour(xs, ys, mo['Z'][0][domain_bounds], levels=np.arange(10, 80, 10), colors='k', zorder=10)
            pylab.contourf(xs, ys, data['t'][wdt][domain_bounds], levels=range(289, 296), cmap=cmap)

            grid.drawPolitical(scale_len=10)

#           pylab.plot(track_xs, track_ys, 'mv-', lw=2.5, mfc='k', ms=8)

            pylab.text(0.05, 0.95, experiments[exp], ha='left', va='top', transform=pylab.gca().transAxes, size=14 * multiplier)
            return
        return doSubplot

    subplots = []
    for exp in experiments.iterkeys():
        subplots.append(subplotFactory(exp, t_ens))

    pylab.figure(figsize=(12, 8))
    pylab.subplots_adjust(left=0.025, bottom=0.1, right=0.875, top=0.975, hspace=0.1, wspace=0.1)

    publicationFigure(subplots, (2, 3), corner='ur')

    cax = pylab.axes((0.90, 0.1125, 0.020, 0.825))
    bar = pylab.colorbar(cax=cax)
    bar.set_ticks(range(287, 296))
    bar.ax.text(3.5, 0.5, r"Temperature (K)", rotation=90, transform=bar.ax.transAxes, size='large', va='center')

    pylab.savefig("cold_pool_%d.png" % t_ens)
    return

if __name__ == "__main__":
    main()
