
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

from datetime import datetime, timedelta

from arpsmodelobs import ARPSModelObsFile
from radarobsfile import RadarObsFile
from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from color_tables import NWSRef
from dataload import loadEnsemble
from computeQuantities import toRecArray
from util import publicationFigure

experiments =     [ '1kmf-sndr0h=25km', '1kmf-zs25-no-05XP', '1kmf-zs25-no-mm-05XP', '1kmf-zs25-no-mm', '1kmf-z-no-snd', '1kmf-z-no-v2' ]
min_ens_members = [ 17,                 31,                  36,                     31,                31,              1              ]

def main():
    base_path = "/caps2/tsupinie/"
    exp_names = { '1kmf-sndr0h=25km':"CTRL", '1kmf-zs25-no-05XP':"NO_MWR", '1kmf-zs25-no-mm-05XP':"NO_MWR_MM", '1kmf-zs25-no-mm':"NO_MM", '1kmf-z-no-snd':"NO_SND", '1kmf-z-no-v2':"NO_V2" }
    experiments =     [ '1kmf-sndr0h=25km', '1kmf-zs25-no-05XP', '1kmf-z-no-snd', '1kmf-zs25-no-mm' ] #, '1kmf-zs25-no-mm-05XP', '1kmf-z-no-v2']
    min_ens_members = [ 17,                 31,                  31,              31                ] #, 36,                     1             ]

    bounds_1sthalf = (slice(105, 160), slice(105, 160))
    bounds_2ndhalf = (slice(130, 185), slice(105, 160))
    grid_1 = goshen_1km_grid(bounds=bounds_1sthalf)
    grid_2 = goshen_1km_grid(bounds=bounds_2ndhalf)
    bounds_1sthalf = grid_1.getBounds()
    bounds_2ndhalf = grid_2.getBounds()
    xs_1, ys_1 = grid_1.getXY()
    xs_2, ys_2 = grid_2.getXY()

    thin_factor = 2
    thin = tuple([slice(None, None, thin_factor)] * 2)

    temp = goshen_1km_temporal(start=14400)
    wind = {}
    for exp, min_ens in zip(experiments, min_ens_members):
        wind[exp] = loadEnsemble("%s%s" % (base_path, exp), [ min_ens ], temp.getTimes(), (['u', 'v', 'w'], toRecArray), {'z':1000}, agl=True)[0]

    def modelSubplotFactory(exp, min_ens, time_sec):
        wdt = temp.getTimes().index(time_sec)
        def doSubplot(multiplier=1.0, layout=(-1, -1)):
            if time_sec < 16200:
                xs, ys = xs_1, ys_1
                domain_bounds = bounds_1sthalf
                grid = grid_1
            else:
                xs, ys = xs_2, ys_2
                domain_bounds = bounds_2ndhalf
                grid = grid_2

            try:
                mo = ARPSModelObsFile("%s/%s/KCYS%03dan%06d" % (base_path, exp, min_ens, time_sec))
            except AssertionError:
                mo = ARPSModelObsFile("%s/%s/KCYS%03dan%06d" % (base_path, exp, min_ens, time_sec), mpi_config=(2, 12))
            except:
                print "Can't load reflectivity ..."
                mo = {'Z':np.zeros((1, 255, 255), dtype=np.float32)}

            pylab.contour(xs, ys, wind[exp]['w'][wdt][domain_bounds], levels=np.arange(2, 102, 2), styles='-', colors='k')
            pylab.contour(xs, ys, wind[exp]['w'][wdt][domain_bounds], levels=np.arange(-100, 0, 2), styles='--', colors='k')

            pylab.quiver(xs[thin], ys[thin], wind[exp]['u'][wdt][domain_bounds][thin], wind[exp]['v'][wdt][domain_bounds][thin])

            pylab.contourf(xs, ys, mo['Z'][0][domain_bounds], levels=np.arange(10, 85, 5), cmap=NWSRef, zorder=-10)

            grid.drawPolitical(scale_len=10)

            row, col = layout
            if col == 1:
                pylab.text(-0.075, 0.5, exp_names[exp], transform=pylab.gca().transAxes, rotation=90, ha='center', va='center', size=12 * multiplier)

        return doSubplot

    def obsSubplotFactory(time):
       def doSubplot(multiplier=1.0, layout=(-1, -1)):
            if (time - datetime(2009, 6, 5, 18, 0, 0)).total_seconds() < 16200:
                xs, ys = xs_1, ys_1
                domain_bounds = bounds_1sthalf
                grid = grid_1
            else:
                xs, ys = xs_2, ys_2
                domain_bounds = bounds_2ndhalf
                grid = grid_2

            try:
                erf = RadarObsFile("qc/1km/KCYS.20090605.%s" % time.strftime("%H%M%S"))
            except:
                print "Can't load reflectivity ..."
                erf = {'Z':np.zeros((1, 255, 255), dtype=np.float32)}

            pylab.contourf(xs, ys, erf['Z'][0][domain_bounds], levels=np.arange(10, 85, 5), cmap=NWSRef)
            grid.drawPolitical(scale_len=10)

            row, col = layout
            if col == 1:
                pylab.text(-0.075, 0.5, "Observations", transform=pylab.gca().transAxes, rotation=90, ha='center', va='center', size=12 * multiplier)

            pylab.text(0.5, 1.075, "%s UTC" % time.strftime("%H%M"), transform=pylab.gca().transAxes, ha='center', va='center', size=12 * multiplier)

       return doSubplot

    pylab.figure(figsize=(18, 21))
    pylab.subplots_adjust(left=0.025, bottom=0.1, right=0.875, top=0.975, hspace=0.05, wspace=0.05)

    subplots = []
    for dt in temp.getDatetimes(aslist=True)[::4]:
        subplots.append(obsSubplotFactory(dt))

    for exp, min_ens in zip(experiments, min_ens_members):
        for time_sec in temp.getTimes()[::4]:
            subplots.append(modelSubplotFactory(exp, min_ens, time_sec))

    publicationFigure(subplots, (5, 4), corner='ur', colorbar=("Reflectivity (dBZ)", "%d", np.arange(10, 85, 5)))

    pylab.savefig("closest/bref.png")
    pylab.close()
    return

if __name__ == "__main__":
    main()
