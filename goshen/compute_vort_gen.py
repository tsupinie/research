
import argparse
import cPickle

import numpy as np
from numpy.lib.recfunctions import append_fields

from computeQuantities import compute3DVorticity, computeVGTensor
from dataload import loadEnsemble, getAxes
from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from interp import getInterpFunctions

def getVGTandUV(**kwargs):
    vg_tensor = computeVGTensor(**kwargs)
    vg_tensor = append_fields(vg_tensor, ['u', 'v'], [ kwargs['u'], kwargs['v'] ])
    return vg_tensor.reshape(kwargs['u'].shape)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp-name', dest='exp_name', required=True)

    args = ap.parse_args()

    np.seterr(all='ignore')
    exp_name = args.exp_name
    base_path = "/caps2/tsupinie/"
    data_path = "%s/%s/" % (base_path, exp_name)

    n_ens_members = 40
    temp = goshen_1km_temporal(start=14400)
    grid = goshen_1km_grid()

    ens_vg = loadEnsemble(data_path, n_ens_members, temp.getTimes(), (['u', 'v', 'w', 'x', 'y', 'z'], getVGTandUV), {'z':500}, agl=True, buffer=True)

    ens_vort = compute3DVorticity(vg_tensor=ens_vg)
    z_stretching = ens_vort['zvort'] * ens_vg['dwdz']
    y_stretching = ens_vort['yvort'] * ens_vg['dvdy']
    x_stretching = ens_vort['xvort'] * ens_vg['dudx']

    horiz_tilting = ens_vort['yvort'] * ens_vg['dwdy'] + ens_vort['xvort'] * ens_vg['dwdx']
    horiz_stretching = (ens_vg['u'] * x_stretching + ens_vg['v'] * y_stretching) / np.hypot(ens_vg['u'], ens_vg['v'])

    cPickle.dump((z_stretching, horiz_tilting, horiz_stretching, ens_vg['u'].mean(axis=0), ens_vg['v'].mean(axis=0)), open("vort_gen_mean_%s.pkl" % exp_name, 'w'), -1)
    return

if __name__ == "__main__":
    main()
