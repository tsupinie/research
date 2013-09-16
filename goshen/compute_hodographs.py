
import numpy as np

from dataload import loadEnsemble
from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from computeQuantities import vectorTuple

import argparse
import cPickle

def main():
    base_path = "/caps2/tsupinie/"
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp-name', dest='exp_name', required=True)

    args = ap.parse_args()

    np.seterr(all='ignore')
    exp_name = args.exp_name
    base_path = "/caps2/tsupinie/"
    data_path = "%s/%s/" % (base_path, exp_name)

    n_ens_members = 40
    temp = goshen_1km_temporal(start=14400, end=18000)
    grid = goshen_1km_grid()

#   Sounding loc: [164124.50758138258, 92544.16037613325]
    y_snd, x_snd = grid.getXY(164, 103)
    print x_snd, y_snd

    ens_hodo = loadEnsemble(data_path, n_ens_members, temp.getTimes(), (['u', 'v'], vectorTuple), {'x':x_snd, 'y':y_snd})

    cPickle.dump(ens_hodo, open("hodo_pkl/%s_hodo.pkl" % exp_name, 'w'), -1)
    return

if __name__ == "__main__":
    main()
