
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from arpsmodelobs import ARPSModelObsFile
from util import publicationFigure

import cPickle
from collections import OrderedDict

def main():
    base_path="/caps2/tsupinie/"
    experiments = OrderedDict([('1kmf-sndr0h=25km', 'CTRL'), ('1kmf-zs25-no-05XP', 'NO_MWR'), ('1kmf-z-no-snd', 'NO_SND'), ('1kmf-zs25-no-mm', 'NO_MM'), ('1kmf-zs25-no-mm-05XP', 'NO_MWR_MM'), ('1kmf-z-no-v2', 'NO_V2')])
    pkl_path = "vort_gen/500m"

    vector_thin = 1
    thin = tuple([ slice(None, None, vector_thin) ] * 2)
    domain_bounds = (slice(110, 135), slice(118, 143))
    grid = goshen_1km_grid(bounds=domain_bounds)
    domain_bounds = grid.getBounds()
    temp = goshen_1km_temporal(start=14400)

    exp_vort = []
    exp_refl = []
    for exp in experiments.iterkeys():
        exp_vort.append(cPickle.load(open("%s/vort_gen_mean_%s.pkl" % (pkl_path, exp), 'r')))
        exp_vort[-1] = exp_vort[-1] + (cPickle.load(open("%s/vort_gen_baroc_mean_%s.pkl" % (pkl_path, exp), 'r')),)

        refl = []
        for time_sec in temp:
            try:
                mo = ARPSModelObsFile("%s/%s/KCYSan%06d" % (base_path, exp, time_sec))
            except AssertionError:
                mo = ARPSModelObsFile("%s/%s/KCYSan%06d" % (base_path, exp, time_sec), mpi_config=(2, 12))
            except:
                mo = {'Z':np.zeros((1, 255, 255), dtype=np.float32)}

            refl.append(mo['Z'][0])
        exp_refl.append(np.array(refl))

    v_stretching, tilting, h_stretching, u_mean, v_mean, baroclinic = zip(*exp_vort)

    v_stretching_max = [ vs.max(axis=0) for vs in v_stretching ]
    tilting_max      = [ tl.max(axis=0) for tl in tilting ]
    h_stretching_max = [ hs.max(axis=0) for hs in h_stretching ]
    baroclinic_max   = [ bc.max(axis=0) for bc in baroclinic ]


    v_stretching = [ vs.mean(axis=0) for vs in v_stretching ]
    tilting      = [ tl.mean(axis=0) for tl in tilting ]
    h_stretching = [ hs.mean(axis=0) for hs in h_stretching ]
    baroclinic   = [ bc.mean(axis=0) for bc in baroclinic ]

    xs, ys = grid.getXY()

    def subplotFactory(exp, vgen, vgen_max, refl, u, v, levels):
        def doSubplot(multiplier=1.0, layout=(-1, -1)):
            pylab.contour(xs, ys, refl[domain_bounds], levels=np.arange(20, 80, 20), colors='#666666')
            pylab.quiver(xs[thin], ys[thin], u[domain_bounds][thin], v[domain_bounds][thin])

            c_level = 4 * len(levels) / 5
            print "Contours at %e s^-2" % (levels[c_level])

            pylab.contour(xs, ys, vgen_max[domain_bounds], levels=[ levels[c_level] ], colors='k', linestyles='-', linewidths=1.5)
            pylab.contour(xs, ys, vgen[domain_bounds], levels=[ levels[c_level] ], colors='k', linestyles='--', linewidths=1.5)

            pylab.contourf(xs, ys, vgen[domain_bounds], levels=levels, cmap=matplotlib.cm.get_cmap('RdBu_r'), zorder=-10)
            grid.drawPolitical()

            pylab.text(0.05, 0.95, experiments[exp], transform=pylab.gca().transAxes, ha='left', va='top', size=14 * multiplier)

            return

        return doSubplot

    levels = np.arange(-0.0001, 0.000105, 0.000005)
    for wdt, time_sec in enumerate(temp):
        subplots = []
        pylab.figure(figsize=(12, 8))
        pylab.subplots_adjust(left=0.05, bottom=0.1, right=0.875, top=0.975, hspace=0.1, wspace=0.1)

        for exp, strt, strt_max, refl, u, v in zip(experiments.iterkeys(), v_stretching, v_stretching_max, exp_refl, u_mean, v_mean):
            subplots.append(subplotFactory(exp, strt[wdt], strt_max[wdt], refl[wdt], u[wdt], v[wdt], levels))

        publicationFigure(subplots, (2, 3), corner='ur', colorbar=(r"Vertical Stretching ($\times$ 10$^3$ s$^{-2}$)", "%.2f", levels[::2], levels[::2] * 1000))

        pylab.savefig("%s/stretching_%d.png" % (pkl_path, time_sec))
        pylab.close()

    levels = np.arange(-0.00007, 0.000075, 0.000005)
    for wdt, time_sec in enumerate(temp):
        subplots = []
        pylab.figure(figsize=(12, 8))
        pylab.subplots_adjust(left=0.05, bottom=0.1, right=0.875, top=0.975, hspace=0.1, wspace=0.1)

        for exp, tilt, tilt_max, refl, u, v in zip(experiments.iterkeys(), tilting, tilting_max, exp_refl, u_mean, v_mean):
            subplots.append(subplotFactory(exp, tilt[wdt], tilt_max[wdt], refl[wdt], u[wdt], v[wdt], levels))

        publicationFigure(subplots, (2, 3), corner='ur', colorbar=(r"Tilting ($\times$ 10$^3$ s$^{-2}$)", "%.2f", levels[::2], levels[::2] * 1000))

        pylab.savefig("%s/tilting_%d.png" % (pkl_path, time_sec))
        pylab.close()

    levels = np.arange(-0.0001, 0.000105, 0.000005)
    for wdt, time_sec in enumerate(temp):
        subplots = []
        pylab.figure(figsize=(12, 8))
        pylab.subplots_adjust(left=0.05, bottom=0.1, right=0.875, top=0.975, hspace=0.1, wspace=0.1)

        for exp, strt, strt_max, refl, u, v in zip(experiments.iterkeys(), h_stretching, h_stretching_max, exp_refl, u_mean, v_mean):
            subplots.append(subplotFactory(exp, strt[wdt], strt_max[wdt], refl[wdt], u[wdt], v[wdt], levels))

        publicationFigure(subplots, (2, 3), corner='ur', colorbar=(r"Horizontal Streamwise Stretching ($\times$ 10$^3$ s$^{-2}$)", "%.2f", levels[::2], levels[::2] * 1000))

        pylab.savefig("%s/horiz_stretching_%d.png" % (pkl_path, time_sec))
        pylab.close()

    levels = np.arange(-0.00005, 0.000055, 0.000005)
    for wdt, time_sec in enumerate(temp):
        subplots = []
        pylab.figure(figsize=(12, 8))
        pylab.subplots_adjust(left=0.05, bottom=0.1, right=0.875, top=0.975, hspace=0.1, wspace=0.1)

        for exp, brcl, brcl_max, refl, u, v in zip(experiments.iterkeys(), baroclinic, baroclinic_max, exp_refl, u_mean, v_mean):
            subplots.append(subplotFactory(exp, brcl[wdt], brcl_max[wdt], refl[wdt], u[wdt], v[wdt], levels))

        publicationFigure(subplots, (2, 3), corner='ur', colorbar=(r"Streamwise Baroclinic Generation ($\times$ 10$^3$ s$^{-2}$)", "%.2f", levels[::2], levels[::2] * 1000))

        pylab.savefig("%s/baroclinic_%d.png" % (pkl_path, time_sec))
        pylab.close()

    return

if __name__ == "__main__":
    main()
