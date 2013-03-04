
import numpy as np

import Nio as nio

from grid import goshen_3km_grid
from util import loadAndInterpolateEnsemble, interpolate, coneHeight
from computeQuantities import computeReflectivity

import glob
from multiprocessing import Process, Pipe, Queue
import cPickle
import gc

def computeRadialVelocity(**kwargs):
    rhat_x, rhat_y, rhat_z = kwargs['rhat']
    return rhat_x * kwargs['u'] + rhat_y * kwargs['v'] + rhat_z * kwargs['w']

def getObsVariables(**kwargs):
    obs = np.empty(kwargs['u'].shape, dtype=[ (v, np.float32) for v in ['u', 'v', 'w', 'pt', 'p', 'qr', 'qs', 'qh'] ])

    obs['u'] = kwargs['u']
    obs['v'] = kwargs['v']
    obs['w'] = kwargs['w']
    obs['pt'] = kwargs['pt']
    obs['p'] = kwargs['p']
    obs['qr'] = kwargs['qr']
    obs['qs'] = kwargs['qs']
    obs['qh'] = kwargs['qh']
    return obs

def computeElevations(ranges, altitudes):
    r_e = 6371000.
    dev_angle = ranges / r_e

    num = (altitudes * np.cos(dev_angle) + (np.cos(dev_angle) - 1) * r_e)
    denom = np.sqrt(altitudes ** 2 + 4 * altitudes * np.sin(dev_angle / 2) ** 2 * r_e - 2 * (np.cos(dev_angle) - 1) * r_e ** 2)
    return np.degrees(np.arcsin(num / denom))

def computeRadialUnitVectors(x_coords, y_coords, radar_elev, levels):
    beam_height = np.empty((len(levels), y_coords.shape[0], x_coords.shape[0]))
    for kdz, radar_tilt in enumerate(levels):
        beam_height[kdz] = coneHeight(np.hypot(*np.ix_(y_coords, x_coords)), radar_tilt, radar_elev)
    
    bh_sq = beam_height ** 2
    yc_sq = y_coords[np.newaxis, :, np.newaxis] ** 2
    xc_sq = x_coords[np.newaxis, np.newaxis, :] ** 2
    slant_dist = np.sqrt(bh_sq + yc_sq + xc_sq)

    return x_coords[np.newaxis, np.newaxis, :] / slant_dist, y_coords[np.newaxis, :, np.newaxis] / slant_dist, beam_height / slant_dist

def beamWeightRadar(obs_dict, radar_unit_vectors, elev_angles, z_coords, levels):
    nz, ny, nx = z_coords.shape
    ens_obs = np.empty((len(radar_unit_vectors.keys()), len(levels), ny, nx), dtype=[('Z', np.float32), ('vr', np.float32)])

    depth = (z_coords[1:] - z_coords[:-1]) / 2.
    depth = np.append(depth, depth[-1].reshape(1, ny, nx), axis=0)

    ens_refl = computeReflectivity(**obs_dict)
    filter_width = 4 * np.log(4.) / 1. ** 2
    for lde, radar in enumerate(sorted(radar_unit_vectors.keys())):
#       radar_elev, radar_x, radar_y = radar_loc[radar]
#       interp_vars = dict([ (f, np.nan * np.ones((len(levels), ny, nx))) for f in obs_vars.dtype.fields.iterkeys() ])

        ens_vel = computeRadialVelocity(rhat=radar_unit_vectors[radar], **obs_dict)

        for index in np.ndindex(ens_obs.shape[1:]):
            radar_tilt = levels[index[0]]

            col_idx = [ slice(None) ]
            col_idx.extend(index[1:])
            col_idx = tuple(col_idx)

            ens_obs['Z'][lde][index] = np.nan
            ens_obs['vr'][lde][index] = np.nan

            sum_idxs = np.where(np.abs(radar_tilt - elev_angles[radar][col_idx]) < 1.0)[0]
            if len(sum_idxs) > 0:
                sum_min, sum_max = sum_idxs.min(), sum_idxs.max()

                sum_idx = [ slice(sum_min, sum_max) ]
                sum_idx.extend(index[1:])
                sum_idx = tuple(sum_idx)

                weight = depth[sum_idx] * np.exp(-filter_width * (radar_tilt - elev_angles[radar][sum_idx]) ** 2)

                ens_obs['Z'][lde][index] = (weight * ens_refl[sum_idx]).sum() / weight.sum()
                ens_obs['vr'][lde][index] = (weight * ens_vel[sum_idx]).sum() / weight.sum()

#           for field in obs_vars.dtype.fields.iterkeys():
#               interp_vars[field][kdz] = interpolate(obs_vars[field][0][0], grid_axes, {'z_base':radar_elev, 'y_base':radar_y, 'x_base':radar_x, 'elev_angle':radar_tilt}, wrap=False)
#               interp_vars[field][kdz] = interpolate(obs_vars[field][0][0], grid_axes, {'z':beam_height[kdz]}, wrap=False)

#       ens_obs['Z'][lde] = computeReflectivity(**interp_vars)
#       ens_obs['vr'][lde] = computeRadialVelocity(rhat=radar_unit_vectors[radar], **interp_vars)
    return ens_obs

def doObsOperator(base_path, grdbas_filename, z_coords, radar_unit_vectors, elev_angles, levels, ens_member, time, comm_pipe=None):
    files = glob.glob("%s/ena%03d.hdf%06d" % (base_path, ens_member, time))

    obs_vars, ens_members, ens_times = loadAndInterpolateEnsemble(files, ['u', 'v', 'w', 'pt', 'p', 'qr', 'qs', 'qh'], getObsVariables, grdbas_filename, prog=False)
    n_ens_members, n_times, nz, ny, nx = obs_vars.shape

    obs_dict = dict([ (k, obs_vars[k][0][0]) for k in obs_vars.dtype.fields.iterkeys() ])

    ens_obs = beamWeightRadar(obs_dict, radar_unit_vectors, elev_angles, z_coords, levels)

    if comm_pipe is not None:
        comm_pipe.put(obs_dict, False)

    cPickle.dump(ens_obs, open("%s/eno%03d.pkl%06d" % (base_path, ens_member, time), 'w'), -1)
    return

def main():
    max_procs = 10
    procs = []
    pipes = {}
    compute_ens_mean = True

    n_ens_members = 40
    t_ens_start = 14400
    t_ens_end = 18000
    dt_ens_step = 300

    radar_loc = {
        'KCYS':(1867, 41.151944, -104.806111),
        'KFTG':(1705, 39.78667, -104.54583),
        'KRIW':(1712, 43.06611, -108.47722),
    }

    grid_obj = goshen_3km_grid()
    base_path = "/caps1/tsupinie/3km-fixed-radar/"
    grdbas_filename = "%s/ena001.hdfgrdbas" % base_path
    levels_vcp11 = [0.5, 0.9, 1.3, 1.8, 2.4, 3.1, 4.0, 5.1, 6.4, 8.0, 10.0, 12.5, 15.6, 19.5] 

    grdbas = nio.open_file(grdbas_filename, mode='r', format='hdf')
    z_coords = grdbas.variables['zp'][:]
    y_coords = grdbas.variables['y'][:]
    x_coords = grdbas.variables['x'][:]
    grid_axes = {'z':z_coords, 'y':y_coords, 'x':x_coords}

    radar_unit_vectors = {}
    radar_loc_xy = {}
    radar_elev_angles = {}

#   beam_height = cPickle.load(open('beam_height.pkl', 'r'))

    for radar, loc in radar_loc.iteritems():
        radar_elev, radar_lat, radar_lon = loc
        radar_x, radar_y = grid_obj(radar_lon, radar_lat)

        radar_loc_xy[radar] = (radar_elev, radar_x, radar_y) 
#       radar_unit_vectors[radar] = computeRadialUnitVectors(x_coords - radar_x, y_coords - radar_y, radar_elev, levels_vcp11)

        ranges = np.hypot(*np.ix_(y_coords - radar_y, x_coords - radar_x))

        bh_sq = z_coords ** 2
        yc_sq = y_coords[np.newaxis, :, np.newaxis] ** 2
        xc_sq = x_coords[np.newaxis, np.newaxis, :] ** 2
        slant_dist = np.sqrt(bh_sq + yc_sq + xc_sq)

        radar_unit_vectors[radar] = (x_coords[np.newaxis, np.newaxis, :] / slant_dist, y_coords[np.newaxis, :, np.newaxis] / slant_dist, z_coords / slant_dist)
        radar_elev_angles[radar] = computeElevations(ranges, z_coords)

    for t_ens in range(t_ens_start, t_ens_end + dt_ens_step, dt_ens_step):
        ensemble_variables = {}

        for n_ens in range(n_ens_members):
            print "Loading ens %03d, time %06d" % (n_ens + 1, t_ens)

            if compute_ens_mean:
                pipe = Queue(len(radar_loc.keys()))
                pipes[n_ens] = pipe
                kwargs = {'comm_pipe':pipe}
            else:
                kwargs = {'comm_pipe':None}

            proc = Process(target=doObsOperator, name="ens=%03d, t=%06d" % (n_ens + 1, t_ens), args=(base_path, grdbas_filename, grid_axes['z'], radar_unit_vectors, radar_elev_angles, levels_vcp11, n_ens + 1, t_ens), kwargs=kwargs)
            proc.start()
            procs.append(proc)

            if len(procs) == max_procs or n_ens + 1 == n_ens_members:
                if compute_ens_mean:
                    for n_ens in sorted(pipes.keys()):
                        obs_vars = pipes[n_ens].get()
                        for var in obs_vars.iterkeys():
                            try:
                                ensemble_variables[var] = (n_ens * ensemble_variables[var] + obs_vars[var]) / (n_ens + 1)
                            except KeyError:
                                ensemble_variables[var] = obs_vars[var]

                    procs = []
                    pipes = {}
                elif len(procs) == max_procs:
                    for proc in procs:
                        proc.join()

                    procs = []

        if compute_ens_mean:
#           for var in ensemble_variables.iterkeys():
#               ensemble_variables[var] = ensemble_variables[var].mean(axis=0)

            ens_mean_obs = beamWeightRadar(ensemble_variables, radar_unit_vectors, radar_elev_angles, grid_axes['z'], levels_vcp11)

            cPickle.dump(ens_mean_obs, open("%s/eomean.pkl%06d" % (base_path, t_ens), 'w'), -1)

        gc.collect()

    for proc in procs:
        proc.join()
    return

if __name__ == "__main__":
    main()
