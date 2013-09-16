
import cPickle
from datetime import datetime, timedelta

import numpy as np
import matplotlib
matplotlib.use('agg')
import pylab
import matplotlib.transforms as transforms

from mpl_toolkits.basemap import Basemap

from legacy import setupMapProjection, goshen_1km_proj, goshen_1km_gs
from util import loadObs

def main():
    exp_names = [ "no-mm", "mm", "mod-05XP" ]
    labels = { "no-mm":"No MM", "mm":"MM", "mod-05XP":"MM + MWR05XP" }
    parameters = [ 't', 'td', 'u', 'v' ]
    param_names = { 't':"Temperature", 'td':"Dewpoint", 'u':r'$u$ Wind', 'v':r'$v$ Wind' }
    units = { 't':r'$^{\circ}$F', 'td':r'$^{\circ}$F', 'u':r'm s$^{-1}$', 'v':r'm s$^{-1}$' }

    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs)
    map = Basemap(**proj)

    times = np.arange(14700, 18300, 300)
    base_time = datetime(2009, 6, 5, 18, 0, 0)
    dt_times = [ base_time + timedelta(seconds=int(t)) for t in times ]
    epoch = datetime(1970, 1, 1, 0, 0, 0)
    base_epoch = (base_time - epoch).total_seconds()

    obs_file_names = ['psu_straka_mesonet.pkl', 'ttu_sticknet.pkl', 'asos.pkl', 'soundings_clip.pkl']
    all_obs = loadObs(obs_file_names, dt_times, map, (goshen_1km_proj['width'], goshen_1km_proj['height']), sounding_obs=['soundings_clip.pkl'])

    ob_nums = []
    for t in times:
        obs_idxs = np.where(all_obs['time'] == base_epoch + t)[0]
        ob_nums.append(len(obs_idxs))

    figures = dict([ (p, pylab.figure()) for p in parameters ])
    axes = {}
    for p in parameters:
        pylab.figure(figures[p].number)
        axes[p] = pylab.axes((0.09, 0.12, 0.82, 0.8))

    for name in exp_names:
        crps = cPickle.load(open("%s_crps.pkl" % name, 'r'))
        for param in parameters:
            pylab.figure(figures[param].number)
            pylab.sca(axes[param])

            pylab.plot(times, crps[param], label=labels[name])

    for param in parameters:
        pylab.figure(figures[param].number)

        num_label_trans = transforms.blended_transform_factory(pylab.gca().transData, pylab.gca().transAxes)

        for t, n_obs in zip(times, ob_nums):
            pylab.text(t, 0.025, "%d" % n_obs, weight='bold', style='italic', size='xx-large', transform=num_label_trans, ha='center', bbox={'facecolor':'#ffffff', 'alpha':0.7})

        pylab.xlim(times.min(), times.max())
        lb_y, ub_y = pylab.ylim()
        pylab.ylim(0, ub_y)

        pylab.xlabel("Time (UTC)", size='large')
        pylab.ylabel("CRPS (%s)" % units[param], size='large')

        pylab.xticks(times, [ (base_time + timedelta(seconds=int(t))).strftime("%H%M") for t in times], rotation=30., size='large')
        pylab.yticks(size='large')

        pylab.legend(loc=1)
        pylab.suptitle("CRPS for %s" % param_names[param])
        pylab.savefig("all_crps_%s.png" % param)
        pylab.close()
    return

if __name__ == "__main__":
    main()
