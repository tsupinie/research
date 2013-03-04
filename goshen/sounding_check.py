
from util import loadObs, loadAndInterpolateEnsemble, setupMapProjection, goshen_1km_proj, goshen_1km_gs
from computeQuantities import theta2Temperature, qv2Dewpoint
from plot_sounding import plotSounding

import numpy as np
import pylab
from mpl_toolkits.basemap import Basemap

from datetime import datetime, timedelta
import glob
from math import floor

def getObsData(**kwargs): 
    obs = np.empty(kwargs['pt'].shape, dtype=[('u', np.float32), ('v', np.float32), ('pt', np.float32), ('p', np.float32), ('qv', np.float32)])
    obs['u'] = kwargs['u']
    obs['v'] = kwargs['v']
    obs['pt'] = kwargs['pt']
    obs['p'] = kwargs['p']
    obs['qv'] = kwargs['qv']

    return obs

def main():
    base_time = datetime(2009, 6, 5, 18, 0, 0)
    epoch = datetime(1970, 1, 1, 0, 0, 0)
    times_seconds = range(14700, 18300, 300)
    times = [ base_time + timedelta(seconds=t) for t in times_seconds ]

    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs)
    map = Basemap(**proj)

    sounding_obs = loadObs(['soundings.pkl'], times, map, sounding_obs=['soundings.pkl'])

    obs_x, obs_y = map(sounding_obs['longitude'], sounding_obs['latitude'])
    obs_z = sounding_obs['elevation']

    start_time = floor(sounding_obs['time'].min() / 300) * 300 - (base_time - epoch).total_seconds()
    sonde_ids = np.unique1d(sounding_obs['id'])

    sondes = {}
    for id in sonde_ids:
        sondes[id] = {'obs':[], 'ens':[] }

    for time in times_seconds[times_seconds.index(start_time):]:
        time_epoch = time + (base_time - epoch).total_seconds()
#       time_base = (epoch + timedelta(seconds=time) - base_time).total_seconds()
        files = glob.glob("/caps1/tsupinie/1km-control-20120712/ena???.hdf%06d" % time)

        round_times = np.round(sounding_obs['time'] / 300) * 300
        time_idxs = np.where(round_times == time_epoch)

        ens_obs, ens_members, ens_times = loadAndInterpolateEnsemble(files, ['u', 'v', 'pt', 'p', 'qv'], getObsData, "/caps1/tsupinie/1km-control-20120712/ena001.hdfgrdbas", 
            {'z':obs_z[time_idxs], 'y':obs_y[time_idxs], 'x':obs_x[time_idxs]}, agl=False, wrap=True)
        
        ens_obs = np.transpose(ens_obs, axes=(2, 0, 1))

        for sonde in sonde_ids:
            sonde_idxs = np.where(sounding_obs['id'][time_idxs] == sonde)

            sondes[sonde]['obs'].extend(sounding_obs[time_idxs[0][sonde_idxs]])
            sondes[sonde]['ens'].extend([ e[:,0] for e in ens_obs[sonde_idxs] ])
    
    for sonde in sonde_ids:
        ens_obs = np.array(sondes[sonde]['ens'], dtype=sondes[sonde]['ens'][0].dtype)
        ens_temp = theta2Temperature(pt=ens_obs['pt'], p=ens_obs['p'])
        ens_dewp = qv2Dewpoint(qv=ens_obs['qv'], p=ens_obs['p'])

        data_obs = np.array(sondes[sonde]['obs'], dtype=sondes[sonde]['obs'][0].dtype)
        order = np.argsort(data_obs['time'])

        time = data_obs['time'][order] - (base_time - epoch).total_seconds()
        obs_temp = data_obs['temp'][order] + 273.15
        obs_dewp = data_obs['dewp'][order] + 273.15

#       pylab.figure(figsize=(8, 10), dpi=100)
#       pylab.axes((0, 0, 1, 1))

        pylab.figure()

        for ens in xrange(ens_obs.shape[1]):
#           plotSounding(None, t=ens_temp[:, ens][order], td=ens_dewp[:, ens][order], p=ens_obs['p'][:, ens][order] / 100., u=ens_obs['u'][:, ens][order], v=ens_obs['v'][:, ens][order])
            pylab.subplot(211)
            pylab.plot(time, ens_temp[:, ens][order], 'r-', linewidth=0.5)
            pylab.plot(time, ens_dewp[:, ens][order], 'g-', linewidth=0.5)

            pylab.subplot(212)
            pylab.plot(time, ens_obs['p'][:, ens][order] / 100., 'b-', linewidth=0.5)

#       plotSounding(None, t=obs_temp, td=obs_dewp, p=data_obs['pres'][order], u=np.ones(order.shape), v=np.zeros(order.shape))
        pylab.subplot(211)
        pylab.plot(time, obs_temp, 'k-', linewidth=1.0)
        pylab.plot(time, obs_dewp, 'k-', linewidth=1.0)

        pylab.subplot(212)
        pylab.plot(time, data_obs['pres'][order], 'k-', linewidth=1.0)

        sonde_name = sonde.replace('/', '_')
        pylab.savefig("sonde_swath_%s.png" % sonde_name)

        pylab.close()       
    return

if __name__ == "__main__":
    main()
