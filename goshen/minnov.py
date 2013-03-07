
import numpy as np

from datetime import datetime, timedelta
import cPickle

from util import runConcurrently
from arpsmodelobs import ARPSModelObsFile
from radarobsfile import RadarObsFile

def minnov(forecast, observation):
    return (forecast - observation).mean()

def doMInnov(radar, model_path, obs_path, t_ens, base_time, analysis=True):
    if analysis:
        state = "an"
    else:
        state = "bg"

    try:
        model_obs = ARPSModelObsFile("%s/%s%s%06d" % (model_path, radar, state, t_ens))
    except AssertionError:
        model_obs = ARPSModelObsFile("%s/%s%s%06d" % (model_path, radar, state, t_ens), mpi_config=(2, 12))
    except IOError:
        return np.nan, np.nan

    radar_obs = RadarObsFile("%s/%s.%s" % (obs_path, radar, (base_time + timedelta(seconds=t_ens)).strftime("%Y%m%d.%H%M%S")))
    good_idxs_refl = np.where((radar_obs['Z'] > -90) & (model_obs['Z'] > -90) & ((radar_obs['Z'] > 15) | (model_obs['Z'] > 15)))
    good_idxs_vel = np.where((radar_obs['vr'] > -90) & (model_obs['vr'] > -90) & ((radar_obs['Z'] > 15) | (model_obs['Z'] > 15)))

    minnov_refl = minnov(model_obs['Z'][good_idxs_refl], radar_obs['Z'][good_idxs_refl])
    minnov_vel = minnov(model_obs['vr'][good_idxs_vel], radar_obs['vr'][good_idxs_vel])

    return minnov_refl, minnov_vel

def main():
    model_paths = [ "/caps1/tsupinie/3km-fixed-radar/", "/caps2/tsupinie/3km-control/", "/caps2/tsupinie/3km-n0r=8e5/", "/caps2/tsupinie/3km-7dBZ,5ms/", "/caps1/tsupinie/3kmf-r0h=12km/" ]
    obs_path = "/data6/tsupinie/goshen/qc/3km/"
    t_ens_start = 10800
    t_ens_stop = 18000
    t_ens_step = 300

    base_time = datetime(2009, 6, 5, 18, 0, 0)

    all_minnov_refl = {}
    all_minnov_vel = {}

    for model_path in model_paths:
        exp_key = model_path.split("/")[-2]
        all_minnov_refl[exp_key] = {}
        all_minnov_vel[exp_key] = {}
        for radar in [ 'KCYS', 'KFTG', 'KRIW' ]:
            minnov_fcst_refl, minnov_fcst_vel = runConcurrently(doMInnov, xrange(t_ens_start, t_ens_stop + t_ens_step, t_ens_step), args=(radar, model_path, obs_path, "__placeholder__", base_time), kwargs={'analysis':False})
            minnov_anal_refl, minnov_anal_vel = runConcurrently(doMInnov, xrange(t_ens_start, t_ens_stop + t_ens_step, t_ens_step), args=(radar, model_path, obs_path, "__placeholder__", base_time), kwargs={'analysis':True})

#           if radar == "KCYS" and model_path == model_paths[-1]:
#               print "fcst refl:", minnov_fcst_refl
#               print "anal refl:", minnov_anal_refl
#               print "times:", times

            all_minnov_refl[exp_key][radar] = np.vstack((minnov_fcst_refl, minnov_anal_refl)).flatten('F')
            all_minnov_vel[exp_key][radar] = np.vstack((minnov_fcst_vel, minnov_anal_vel)).flatten('F')

    cPickle.dump(all_minnov_refl, open("all_minnov_refl.pkl", 'w'), -1)
    cPickle.dump(all_minnov_vel, open("all_minnov_vel.pkl", 'w'), -1)

    return

if __name__ == "__main__":
   main()
