
import numpy as np

from dataload import loadEnsemble
from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from computeQuantities import toRecArray

import argparse
import cPickle

def createDeviations(args):
    np.seterr(all='ignore')
    exp_name = args.exp_name
    base_path = "/caps2/tsupinie/"
    data_path = "%s/%s/" % (base_path, exp_name)

    n_ens_members = 40
    temp = goshen_1km_temporal(start=args.ens_time, end=args.ens_time)
    grid = goshen_1km_grid()

    ens = loadEnsemble(data_path, n_ens_members, temp.getTimes(), ([ 'u', 'v', 'pt' ], toRecArray))

    ens_dev = np.empty(ens.shape[:1], dtype=ens.dtype)

    for field in ens.dtype.fields.iterkeys():
        ens_mean = ens[field].mean(axis=0)
        print field
        for lde in range(ens.shape[0]):
            mem_dev = np.sqrt(np.mean((ens[field][lde] - ens_mean) ** 2)) / np.std(ens_mean, ddof=1)

    cPickle.dump(ens_dev, open("closest/ens_dev_%s.pkl" % exp_name, 'w'), -1)
    return

def compileDeviations():
    experiments = ['1kmf-sndr0h=25km', '1kmf-zs25-no-05XP', '1kmf-zs25-no-mm-05XP', '1kmf-zs25-no-mm', '1kmf-z-no-snd', '1kmf-z-no-v2']
    temp = goshen_1km_temporal(start=14400, end=14400)

    deviations = []
    for exp in experiments:
        exp_dev = []
        for t_ens in temp:
            exp_dev.append(cPickle.load(open("closest/ens_dev_%s_%d.pkl" % (exp, t_ens), 'r')))

        deviations.append(np.array(exp_dev, dtype=exp_dev[0].dtype))
    deviations = np.array(deviations, dtype=deviations[0].dtype)
    print deviations.shape
    deviations = np.array([ sum(deviations[idx]) / len(deviations[idx]) for idx in np.ndindex(deviations.shape) ]).reshape(deviations.shape).mean(axis=1)
    min_ens_member = np.argmin(deviations, axis=1)
    for exp, min_ens in zip(experiments, min_ens_member + 1):
        print "Closest ensemble member for %s: %d" % (exp, min_ens)
    return

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp-name', dest='exp_name')
    ap.add_argument('--ens-time', dest='ens_time', type=int)

    args = ap.parse_args()
    if args.exp_name is not None and args.ens_time is not None:
        createDeviations(args)
    else:
        compileDeviations()

    return

if __name__ == "__main__":
    main()
