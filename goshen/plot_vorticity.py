
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from util import publicationFigure
from arpsmodelobs import ARPSModelObsFile

import cPickle
from collections import OrderedDict

def main():
    base_path = "/caps2/tsupinie/"
#   experiments = OrderedDict([('1kmf-sndr0h=25km', 'CTRL'), ('1kmf-zs25-no-05XP', 'NO_MWR'), ('1kmf-z-no-snd', 'NO_SND'), ('1kmf-zs25-no-mm', 'NO_MM'), ('1kmf-zs25-no-mm-05XP', 'NO_MWR_MM'), ('1kmf-z-no-v2', 'NO_V2')])
    experiments = OrderedDict([('1kmf-z04vr=30dBZ', 'MWR_VR_THRESH'), ('1kmf-sndr0h=25km', 'CTRL'), ('1kmf-zs25-no-05XP', 'NO_MWR')])

    domain_bounds = (slice(110, 135), slice(118, 143))
    grid = goshen_1km_grid(bounds=domain_bounds)
    domain_bounds = grid.getBounds()
    temp = goshen_1km_temporal(start=14400, end=14400)

    xs, ys = grid.getXY()
    levels = np.arange(-0.030, 0.033, 0.003)

    exp_vort = []
    exp_refl = []

    for exp in experiments.iterkeys():
        vort = cPickle.load(open("vort_pkl/vorticity_%s.pkl" % exp, 'r'))

        refl = []
        for time_sec in temp:        
            try:
                mo = ARPSModelObsFile("%s/%s/KCYSan%06d" % (base_path, exp, time_sec))
            except AssertionError:
                mo = ARPSModelObsFile("%s/%s/KCYSan%06d" % (base_path, exp, time_sec), mpi_config=(2, 12))
            except:
                mo = {'Z':np.zeros((1, 255, 255), dtype=np.float32)}
            refl.append(mo)

        if exp not in [ "1kmf-zs25-no-KCYSvr", "1kmf-z04vr=30dBZ" ]:
            exp_vort.append(vort[:, 4:5])
        else:
            exp_vort.append(vort)

        exp_refl.append(np.array(refl))

    def subplotFactory(exp, exp_vort, exp_refl):
        def doSubplot(multiplier=1.0, layout=(-1, -1)):

            pylab.quiver(xs, ys, exp_vort['u'].mean(axis=0)[domain_bounds], exp_vort['v'].mean(axis=0)[domain_bounds])
            pylab.contour(xs, ys, exp_refl['Z'][0][domain_bounds], colors='#666666', levels=np.arange(20, 80, 20))

            pylab.contour(xs, ys, exp_vort['vort'].mean(axis=0)[domain_bounds], colors='k', linestyles='--', linewidths=1.5, levels=[ 0.015 ])
            pylab.contour(xs, ys, exp_vort['vort'].max(axis=0)[domain_bounds], colors='k', linestyles='-', linewidths=1.5, levels=[ 0.015 ])

            pylab.contourf(xs, ys, exp_vort['vort'].mean(axis=0)[domain_bounds], cmap=matplotlib.cm.get_cmap('RdBu_r'), levels=levels, zorder=-10)
            grid.drawPolitical()

            pylab.text(0.05, 0.95, experiments[exp], ha='left', va='top', transform=pylab.gca().transAxes, size=14 * multiplier)
            return

        return doSubplot

    for wdt, time_sec in enumerate(temp):
        print time_sec
        subplots = []
        for exp_name, vort, refl in zip(experiments.iterkeys(), exp_vort, exp_refl):
            print vort.shape
            subplots.append(subplotFactory(exp_name, vort[:, wdt], refl[wdt]))

        pylab.figure(figsize=(12, 8))
        pylab.subplots_adjust(left=0.05, bottom=0.1, right=0.875, top=0.975, hspace=0.1, wspace=0.1)
        publicationFigure(subplots, (1, 3), corner='ur', colorbar=(r'Vorticity ($\times$ 10$^3$ s$^{-1}$)', "%d", levels, np.round(1000 * levels)))
        pylab.savefig("vorticity_z04vrcompare_%d.png" % time_sec)
    return

if __name__ == "__main__":
    main()
