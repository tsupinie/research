
from legacy import loadAndInterpolateEnsemble, setupMapProjection, goshen_1km_proj, goshen_1km_gs, drawPolitical
from util import loadObs, probMatchMean
from computeQuantities import theta2Temperature, qv2Dewpoint, computeReflectivity
from plot_sounding import plotSkewTBackground, plotProfile, plotWinds

import pylab
from mpl_toolkits.basemap import Basemap
import numpy as np

import glob
import cPickle
from datetime import datetime, timedelta

def getObsData(**kwargs): 
    obs = np.empty(kwargs['pt'].shape, dtype=[('u', np.float32), ('v', np.float32), ('pt', np.float32), ('p', np.float32), ('qv', np.float32)])
    obs['u'] = kwargs['u']
    obs['v'] = kwargs['v']
    obs['pt'] = kwargs['pt']
    obs['p'] = kwargs['p']
    obs['qv'] = kwargs['qv']

    return obs

def getReleaseCoords(sounding):
    release_loc = sounding['release_loc'].split(", ")
    return float(release_loc[2]), float(release_loc[3])

def fcstSounding(ens_obs, sounding_obs, file_name):
    ens_t = theta2Temperature(pt=ens_obs['pt'], p=ens_obs['p'])
    ens_td = qv2Dewpoint(qv=ens_obs['qv'], p=ens_obs['p'])
    ens_p = ens_obs['p']

    pylab.figure(figsize=(8, 10))
    plotSkewTBackground(pylab.gca())

#   print ens_t.shape

    for lde in xrange(ens_obs.shape[0]):
#       print ens_t[lde] - 273.15
#       print ens_p[lde] / 100.

        order = np.argsort(ens_p[lde])

#       print ens_t[lde][order].shape, ens_p[lde][order].shape, ens_td[lde][order].shape
#       print ens_t[lde][order]
#       print ens_p[lde][order]

        plotProfile(ens_t[lde][order] - 273.15, ens_p[lde][order] / 100., color='r', lw=0.5)
        plotProfile(ens_td[lde][order] - 273.15, ens_p[lde][order] / 100., color='g', lw=0.5)

    order = np.argsort(sounding_obs['pres'])
    temp_C = 5. / 9. * (sounding_obs['temp'][order] - 32)
    dewp_C = 5. / 9. * (sounding_obs['dewp'][order] - 32)

    plotProfile(temp_C, sounding_obs['pres'][order], color='k', lw=1.5)
    plotProfile(dewp_C, sounding_obs['pres'][order], color='k', lw=1.5)

    pylab.savefig(file_name)
    pylab.close()
    return

def fcstProfile(ens_base, sounding, map, file_name):
    u_obs_snd = sounding['u_wind']
    v_obs_snd = sounding['v_wind']
    t_obs_snd = sounding['temperature']
    td_obs_snd = sounding['dewpoint']
    p_obs_snd = sounding['pressure']

    good = (u_obs_snd != 9999.0) & (v_obs_snd != 9999.0) & (t_obs_snd != 999.0) & (td_obs_snd != 999.0) & (p_obs_snd != 9999.0)
    good_idxs = np.where(good)[0]

    release_lon, release_lat = getReleaseCoords(sounding)

    snd_x, snd_y = map(release_lon, release_lat)

#   files = glob.glob("%s/enf???.hdf014400" % ens_base)
#   ens_obs, ens_members, ens_times = loadAndInterpolateEnsemble(files, ['u', 'v', 'pt', 'p', 'qv'], getObsData, "%s/enf001.hdfgrdbas" % ens_base, 
#      {'y':snd_y, 'x':snd_x}, )

#   mean_pt = ens_obs['pt'][:, 0, :].mean(axis=0)
#   mean_p  = ens_obs['p' ][:, 0, :].mean(axis=0)
#   mean_qv = ens_obs['qv'][:, 0, :].mean(axis=0)
#   mean_u  = ens_obs['u' ][:, 0, :].mean(axis=0)
#   mean_v  = ens_obs['v' ][:, 0, :].mean(axis=0)

#   mean_t  = theta2Temperature(pt=mean_pt, p=mean_p) - 273.15
#   mean_td = qv2Dewpoint(qv=mean_qv, p=mean_p) - 273.15

#   plotSounding("fcst.bushnell.014400.png", t=mean_t - 273.15, p=(mean_p / 100.), td=mean_td - 273.15, u=mean_u, v=mean_v)
    pylab.figure(figsize=(8, 10))

    plotSkewTBackground(pylab.gca())

#   ens_p = ens_obs['p'][:, 0, :]
#   ens_t  = theta2Temperature(pt=ens_obs['pt'][:, 0, :], p=ens_p) - 273.15
#   ens_td = qv2Dewpoint(qv=ens_obs['qv'][:, 0, :], p=ens_p) - 273.15

#   for lde, ens in enumerate(ens_members):
#       plotProfile(ens_t[lde], ens_p[lde] / 100., color='r', lw=0.5)
#       plotProfile(ens_td[lde], ens_p[lde] / 100., color='g', lw=0.5)

    plotProfile(t_obs_snd[good_idxs], p_obs_snd[good_idxs], color='r', lw=1.5)
    plotProfile(td_obs_snd[good_idxs], p_obs_snd[good_idxs], color='g', lw=1.5)
    plotWinds(u_obs_snd[good_idxs], v_obs_snd[good_idxs], p_obs_snd[good_idxs], color='k')

    pylab.savefig(file_name)
    pylab.close()
    return

def getObsData(**kwargs): 
    obs = np.empty(kwargs['pt'].shape, dtype=[('u', np.float32), ('v', np.float32), ('pt', np.float32), ('p', np.float32), ('qv', np.float32)])
    obs['u'] = kwargs['u']
    obs['v'] = kwargs['v']
    obs['pt'] = kwargs['pt']
    obs['p'] = kwargs['p']
    obs['qv'] = kwargs['qv']

    return obs

def main():
    exp_name = "mod-05XP"
    base_path = "/caps1/tsupinie/1km-control-%s" % exp_name
    files = glob.glob("%s/ena???.hdf018000" % base_path)
    radar_elev, radar_lat, radar_lon = 1883, 41.151944, -104.806111

    nicknames = ['fcst.wynssl.014400.png', 'fcst.bushnell.014400.png', 'fcst.ncar2.014400.png']

    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs)
    map = Basemap(**proj)

#   radar_x, radar_y = map(radar_lon, radar_lat)

#   coords = []

#   base_time = datetime(2009, 6, 5, 18, 0, 0)
#   soundings = loadObs(['soundings_clip.pkl'], [ base_time + timedelta(seconds=18000) ], map, (goshen_1km_proj['width'], goshen_1km_proj['height']), sounding_obs=['soundings_clip.pkl'])

#   sounding_ids = np.unique1d(soundings['id'])

#   soundings = soundings[np.where(soundings['id'] == "Bush")]

#   obs_x, obs_y = map(soundings['longitude'], soundings['latitude'])    
#   obs_p = soundings['pres'] * 100

#   ens_obs, ens_members, ens_times = loadAndInterpolateEnsemble(files, ['u', 'v', 'pt', 'p', 'qv'], getObsData, "/caps1/tsupinie/1km-control-mod-05XP/ena001.hdfgrdbas", 
#       {'z':obs_p, 'y':obs_y, 'x':obs_x}, agl=False, wrap=True, coords='pres')

#   print ens_obs.shape

#   for snd_id in ['Bush']: # sounding_ids:
#       sounding_idxs = np.where(soundings['id'] == snd_id)[0]

#       snd_nickname = snd_id.replace('/', '_')
#       fcstSounding(ens_obs[:, 0, sounding_idxs], soundings[sounding_idxs], "fcst.%s.png" % snd_nickname)

    for snd_file, plot_file in zip(['proximity1.pkl', 'proximity2.pkl', 'proximity3.pkl'], nicknames):
        obs_sounding = cPickle.load(open(snd_file, 'r'))
        fcstProfile(base_path, obs_sounding, map, plot_file)
#       coords.append(getReleaseCoords(obs_sounding))

### Plot the locations of the sounding obs on the probability-matched reflectivity
#   ens_refl, ens_members, ens_times = loadAndInterpolateEnsemble(files, ['pt', 'p', 'qr', 'qs', 'qh'], computeReflectivity, "%s/ena001.hdfgrdbas" % base_path, 
#       {'z_base':radar_elev, 'y_base':radar_y, 'x_base':radar_x, 'elev_angle':0.5}, agl=False, wrap=True) #, aggregator=lambda x: np.mean(x, axis=0))

#   refl_ens_mean = probMatchMean(ens_refl)

#   print ens_refl.shape
#   print refl_ens_mean.shape

#   pylab.figure()

#   dx, dy = goshen_1km_gs
#   nx, ny = refl_ens_mean[0].shape
#   xs, ys = np.meshgrid(dx * np.arange(nx), dy * np.arange(ny))
#   pylab.contourf(xs, ys, refl_ens_mean[0], levels=np.arange(10, 80, 10))

#   for snd_id, marker in zip(sounding_ids, ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'L', 'M'][:len(sounding_ids)]):
#       sounding_idxs = np.where(soundings['id'] == snd_id)[0]

#       lats = soundings['latitude'][sounding_idxs]
#       lons = soundings['longitude'][sounding_idxs]

#       order = np.argsort(soundings['pres'][sounding_idxs])
#       snd_xs, snd_ys = map(lons[order], lats[order])

#       pylab.plot(snd_xs, snd_ys, 'ko-', markersize=2)
#       if snd_id == "88 a":
#           pylab.text(snd_xs[-1] - 6000, snd_ys[-1] + 2000, marker, fontsize='x-large', fontweight='bold', ha='right', va='bottom', bbox={'facecolor':'w', 'alpha':0.7})
#       else:
#           pylab.text(snd_xs[-1] - 2000, snd_ys[-1] - 2000, marker, fontsize='x-large', fontweight='bold', ha='right', va='top', bbox={'facecolor':'w', 'alpha':0.7})

##   for coord, name in zip(coords, nicknames):
##       snd_x, snd_y = map(*coord)
##       pylab.plot(snd_x, snd_y, 'ko')
##       pylab.text(snd_x + 1000, snd_y + 1000, name, ha='left', va='bottom')

#   drawPolitical(map)

##  pylab.savefig("fcst_snd_locations.png")
#   pylab.savefig("snd_locations.png")
#   pylab.close()
    
    return

if __name__ == "__main__":
    main()
