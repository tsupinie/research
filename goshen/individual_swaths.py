
import numpy as np

import glob
import cPickle

from util import loadAndInterpolateRun
from computeQuantities import computeReflectivity

def main():
    cmix_params = [ 0.0010, 0.0015, 0.0020, 0.0025, 0.0030 ]
    interp_height = 975

    refl_vars = [ 'p', 'pt', 'qr', 'qs', 'qh' ]

    param_runs = {}
    max_runs = {}

    for cmix in cmix_params:
        data_dir = "cmix=%6.4f" % cmix

        files = glob.glob("%s/250m-det.hdf0?????" % data_dir)
        grdbas_file = "%s/250m-det.hdfgrdbas" % data_dir

        param_runs[data_dir], times = loadAndInterpolateRun(files, interp_height, refl_vars, computeReflectivity, grdbas_file=grdbas_file)
        print param_runs[data_dir].shape
        max_runs[data_dir] = param_runs[data_dir].max(axis=0)

    cPickle.dump(max_runs, open("250m-det.pkl", 'w'), -1)
    return

if __name__ == "__main__":
    main()
