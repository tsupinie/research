
import cPickle

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

from util import loadObs
from grid import goshen_1km_grid
from temporal import goshen_1km_temporal

def main():
    experiments = [ '1kmf-sndr0h=50km', '1kmf-zs-no-05XP', '1kmf-zs-no-mm-05XP', '1kmf-zs-no-mm', '1kmf-z-no-snd', '1kmf-z-no-v2' ]

    grid = goshen_1km_grid(bounds=(slice(100, 180), slice(90, 170)))
    temp = goshen_1km_temporal(start=14400)

    obs_file_names = ['psu_straka_mesonet.pkl', 'ttu_sticknet.pkl', 'asos.pkl']
    all_obs = loadObs(obs_file_names, temp.getDatetimes(aslist=True), grid, grid.getWidthHeight())

    ens_obs = {}
    for exp in experiments:
        ens_obs[exp] = cPickle.load(open("cold_pool_obs_%s.pkl" % exp, 'r'))

    mm_ids = np.unique1d(all_obs['id'])

    for id in mm_ids:
        id_idxs = np.where(all_obs['id'] == id)
        for ob_var, ens_var in [('temp', 't'), ('dewp', 'td')]:

            pylab.figure()
            pylab.plot(all_obs['time'][id_idxs], 5 / 9. * (all_obs[ob_var][id_idxs] - 32), 'k-', label='Observed')

            for exp_name, exp in ens_obs.iteritems():
                pylab.plot(all_obs['time'][id_idxs], exp[ens_var][id_idxs] - 273.15, label=exp_name)

            pylab.xticks(temp.getEpochs(aslist=True), temp.getStrings("%H%M", aslist=True), rotation=30)
            pylab.xlim(temp.getEpochs(aslist=True)[0], temp.getEpochs(aslist=True)[-1])
            pylab.legend(loc=1)

            pylab.savefig("mm_timeseries_%s_%s.png" % (ens_var, id))
            pylab.close()        

    return

if __name__ == "__main__":
    main()
