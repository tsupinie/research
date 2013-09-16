
import numpy as np

import Nio as nio

import matplotlib
matplotlib.use('agg')
import pylab
import matplotlib.cm as cm

import cPickle
from datetime import datetime,timedelta
from collections import OrderedDict

from temporal import goshen_1km_temporal
from util import publicationFigure

def findHeights(grdbas, bounds):
    hdf = nio.open_file(grdbas, mode='r', format='hdf')

    bounds_x, bounds_y = bounds
    column_x = (bounds_x.start + bounds_x.stop) / 2
    column_y = (bounds_y.start + bounds_y.stop) / 2

    return hdf.variables['zp'][:, column_y, column_x]

def boundCoordinate(coords):
    d_coord = (coords[1:] - coords[:-1]) / 2
    lb_coord = coords[0] - d_coord[0]
    ub_coord = coords[-1] + d_coord[-1]

    return np.append(np.append([ lb_coord ], coords[1:] - d_coord), ub_coord)

def main():
#   exp_names = [ 'sndr0h=25km', 'zs25-no-05XP', 'z-no-snd', 'zs25-no-mm', 'zs25-no-mm-05XP', 'z-no-v2' ]
#   exp_names = [ 'sndr0h=25km', 'zs25-offtime-05XP' ]
    experiments = OrderedDict([('1kmf-sndr0h=25km', 'CTRL'), ('1kmf-zs25-no-05XP', 'NO_MWR'), ('1kmf-z-no-snd', 'NO_SND'), ('1kmf-zs25-no-mm', 'NO_MM')]) #, ('1kmf-zs25-no-mm-05XP', 'NO_MWR_MM'), ('1kmf-z-no-v2', 'NO_V2')])
    temp = goshen_1km_temporal(start=14400)    

    z_column = findHeights("/caps2/tsupinie/%s/ena001.hdfgrdbas" % experiments.keys()[0], (slice(90, 170), slice(100, 180)))
    cutoff = 14

    vort_data = {}
    max_vorticity = 0.035
    min_vorticity = 0.0075

    def subplotFactory(exp):
        def doSubplot(multiplier=1.0, layout=(-1, -1)):
            exp_vort = cPickle.load(open("vort_time_height_%s.pkl" % exp[5:], 'r'))
            vort_data[exp] = exp_vort
            print "Max vorticity:", exp_vort.max()

            cmap = cm.get_cmap('RdYlBu_r')
            cmap.set_under('#ffffff')

#           pylab.pcolormesh(boundCoordinate(np.array(temp.getTimes())), boundCoordinate(z_column / 1000.), exp_vort.T, cmap=cmap, vmin=min_vorticity, vmax=max_vorticity)
            pylab.contourf(temp.getTimes(), z_column / 1000., exp_vort.T, cmap=cmap, levels=np.arange(min_vorticity, max_vorticity + 0.0024, 0.0025))

            pylab.xlim(temp.getTimes()[0], temp.getTimes()[-1])
            pylab.ylim([z_column[0] / 1000, 10])

            layout_r, layout_c = layout
            if layout_r == 2 and layout_c == 1 or layout_r == -1 and layout_c == -1:
                pylab.xlabel("Time (UTC)", size='large')
                pylab.ylabel("Height (km MSL)", size='large')

#           pylab.axvline(14400, color='k', linestyle=':')

                pylab.xticks(temp.getTimes(), temp.getStrings('%H%M', aslist=True), size='large')
                pylab.yticks(size='large')
            else:
                pylab.xticks(temp.getTimes(), [ '' for t in temp ])
                pylab.yticks(pylab.yticks()[0], [ '' for z in pylab.yticks()[0] ])

            pylab.text(0.05, 0.95, experiments[exp], ha='left', va='top', transform=pylab.gca().transAxes, size=18 * multiplier)
            pylab.xlim(temp.getTimes()[0], temp.getTimes()[3])
            pylab.ylim([z_column[1] / 1000, 6])
#           pylab.title(exp)
            return

        return doSubplot

    def diffSubplotFactory(exp):
        def doSubplot(multiplier=1.0, layout=(-1, -1)):
            exp_vort = vort_data[exp]
            pylab.pcolormesh(boundCoordinate(np.array(temp.getTimes())), boundCoordinate(z_column / 1000.), exp_vort.T - vort_data['sndr0h=50km'].T, 
                cmap=cm.get_cmap('RdBu_r'), vmin=-max_vorticity / 2, vmax=max_vorticity / 2)

            layout_r, layout_c = layout
            if layout_r == 2 and layout_c == 1 or layout_r == -1 and layout_c == -1:
                pylab.xlabel("Time (UTC)", size='large')
                pylab.ylabel("Height (km MSL)", size='large')

            pylab.xticks(temp.getTimes(), temp.getStrings('%H%M', aslist=True), size='large')
            pylab.yticks(size='large')

            pylab.xlim(temp.getTimes()[0], temp.getTimes()[3])
            pylab.ylim([z_column[1] / 1000, 6])

#           pylab.title(exp)
            return

        return doSubplot

    subplots = []
    diff_subplots = []
    for exp in experiments.keys():
        subplots.append(subplotFactory(exp))
#       if exp != 'sndr0h=50km':
#           diff_subplots.append(diffSubplotFactory(exp))

#   pylab.figure(figsize=(12, 9))
    pylab.figure(figsize=(8.5, 9))
    pylab.subplots_adjust(left=0.05, bottom=0.1, right=0.875, top=0.975, hspace=0.1, wspace=0.1)

    ticks = np.arange(min_vorticity, max_vorticity + 0.0024, 0.0025)
    publicationFigure(subplots, (2, 2), corner='ur', colorbar=(r"Vertical Vorticity ($\times$ 10$^{-3}$ s$^{-1}$)", '%.1f', ticks, 1000 * ticks))

#   pylab.suptitle("Ensemble Mean of Max $\zeta$")
    pylab.savefig("vort_time_height_four.png")
    pylab.close()

    return

if __name__ == "__main__":
    main()
