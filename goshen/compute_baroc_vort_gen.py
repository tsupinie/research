
import numpy as np

import cPickle
import argparse

from computeQuantities import computeBaroclinicVortGen
from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from dataload import loadEnsemble

def main():
    np.seterr(all='ignore')

    ap = argparse.ArgumentParser()
    ap.add_argument('--exp-name', dest='exp_name', required=True)

    args = ap.parse_args()

    exp_name = args.exp_name
    base_path = "/caps2/tsupinie/"
    data_path = "%s/%s/" % (base_path, exp_name)

    n_ens_members = 40
    temp = goshen_1km_temporal(start=14400)
    grid = goshen_1km_grid()

    ens_bvg = loadEnsemble(data_path, n_ens_members, temp.getTimes(), (['u', 'v', 'pt', 'p', 'qv', 'qc', 'qi', 'qr', 'qs', 'qh', 'x', 'y'], computeBaroclinicVortGen), {'z':500}, agl=True)

    cPickle.dump(ens_bvg, open("vort_gen_baroc_mean_%s.pkl" % exp_name, 'w'), -1)
    return

if __name__ == "__main__":
    main()
