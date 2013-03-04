
import numpy as np

import Nio as nio

import matplotlib
matplotlib.use('agg')
import pylab

import cPickle
from datetime import datetime,timedelta

def findHeights(grdbas, bounds):
    hdf = nio.open_file(grdbas, mode='r', format='hdf')

    bounds_x, bounds_y = bounds
    column_x = (bounds_x.start + bounds_x.stop) / 2
    column_y = (bounds_y.start + bounds_y.stop) / 2

    return hdf.variables['zp'][:, column_y, column_x]

def main():
    exp_names = [ "mod-05XP", "mm", "no-mm" ]
    times = np.arange(10800, 18300, 300)
    base_time = datetime(2009, 6, 5, 18, 0, 0)

    z_column = findHeights("/caps1/tsupinie/1km-control-mod-05XP/ena001.hdfgrdbas", (slice(90, 170), slice(100, 180)))
    cutoff = 14

    vort_data = {}
    for exp in exp_names:
        exp_vort = cPickle.load(open("vort_time_height_%s_anal.pkl" % exp, 'r'))
        exp_vort = np.concatenate((exp_vort, cPickle.load(open("vort_time_height_%s_fcst.pkl" % exp, 'r'))))

        exp_vort = exp_vort[:cutoff]

        pylab.figure() #figsize=(12, 6))
        pylab.axes((0.1, 0.125, 0.8, 0.8))

        pylab.contourf(times[:cutoff], z_column / 1000., exp_vort[:cutoff].T, cmap=matplotlib.cm.get_cmap('Reds'), levels=np.arange(0, 0.036, 0.003))
        pylab.colorbar()

        pylab.xlabel("Time (UTC)", size='large')
        pylab.ylabel("Height (km)", size='large')

        pylab.axvline(14400, color='k', linestyle=':')

        pylab.xticks(times[:cutoff], [ (base_time + timedelta(seconds=int(t))).strftime("%H%M") for t in times[:cutoff] ], rotation=30, size='large')
        pylab.yticks(size='large')

        y_lb, y_ub = pylab.ylim()
        pylab.ylim([y_lb, 10])

        pylab.suptitle("Time-Height plot of maximum vertical vorticity")
        pylab.savefig("vort_time_height_%s.png" % exp)
        pylab.close()

    return

if __name__ == "__main__":
    main()
