
from util import loadAndInterpolateEnsemble
from computeQuantities import computeVorticity

import matplotlib
matplotlib.use('agg')
import pylab
import numpy as np

import glob
import argparse
import cPickle

def plotTimeHeight(time_height, times, title, file_name):
    pylab.figure()

    levels = time_height.shape[1]
    zs, ts = np.meshgrid(np.arange(levels), times)

    pylab.contourf(ts, zs, time_height, cmap=matplotlib.cm.get_cmap('Reds'), levels=np.arange(0, 0.036, 0.003))
    pylab.colorbar()

    pylab.xlabel("Time (s)")
    pylab.ylabel("Model Vertical Level")

    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp-name', dest='exp_name', required=True)

    args = ap.parse_args()

    bounds = (slice(90, 170), slice(100, 180))
    base_path = "/caps1/tsupinie/1km-control-%s" % args.exp_name

#   files = glob.glob("%s/ena???.hdf014[47]00" % base_path)
#   files.extend(glob.glob("%s/ena???.hdf01[5678]?00" % base_path))

    files = glob.glob("%s/ena???.hdf01[0123]*" % base_path)
    files.extend(glob.glob("%s/ena???.hdf014100" % base_path))

    ens_vort, ens_members, ens_times = loadAndInterpolateEnsemble(files, ['u', 'v', 'dx', 'dy'], computeVorticity, "%s/ena001.hdfgrdbas" % base_path)

    full_bounds = [ slice(None) ] * 3
    full_bounds.extend(bounds)

    print ens_vort.shape
    time_height = ens_vort[full_bounds].mean(axis=0).max(axis=-1).max(axis=-1)
    print time_height.shape

    cPickle.dump(time_height, open("vort_time_height_%s_anal.pkl" % args.exp_name, 'w'), -1)

    plotTimeHeight(time_height, ens_times, "Time-Height plot of maximum vertical vorticity", "vort_time_height_%s.png" % args.exp_name)

    return

if __name__ == "__main__":
    main()
