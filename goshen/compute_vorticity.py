
import numpy as np

from dataload import loadEnsemble
from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from computeQuantities import computeVorticity

import cPickle
import argparse

def getVortUV(**kwargs):
    vortuv = np.empty(kwargs['u'].shape, dtype=[('vort', np.float32), ('u', np.float32), ('v', np.float32)])

    vortuv['vort'] = computeVorticity(**kwargs)
    vortuv['u'] = kwargs['u']
    vortuv['v'] = kwargs['v']

    return vortuv

def main():
    base_path = "/caps2/tsupinie/"
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp-name', dest='exp_name', required=True)
    ap.add_argument('--fcst', dest='fcst', action='store_true')

    args = ap.parse_args()

    np.seterr(all='ignore')
    exp_name = args.exp_name
    base_path = "/caps2/tsupinie/"
    data_path = "%s/%s/" % (base_path, exp_name)

    n_ens_members = 40
    temp = goshen_1km_temporal(start=14400, end=14400)
    grid = goshen_1km_grid()

    ens_vort = loadEnsemble(data_path, n_ens_members, temp.getTimes(), (['u', 'v', 'dx', 'dy'], getVortUV), {'z':1000}, agl=True, fcst=args.fcst)
    
    print ens_vort.shape

    fcst_str=""
    if args.fcst: fcst_str = "_fcst"

    cPickle.dump(ens_vort, open("vort_pkl/vorticity%s_%s.pkl" % (fcst_str, exp_name), 'w'), -1)
    return

if __name__ == "__main__":
    main()
