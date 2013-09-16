
from legacy import loadAndInterpolateEnsemble, setupMapProjection, interpolate, goshen_1km_proj, goshen_1km_gs, _makeZCoordsAGL
from util import decompressVariable, loadObs

import pylab
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import CirclePolygon, Rectangle, PathPatch, Polygon
from matplotlib.path import Path

import numpy as np
import Nio as nio
from scipy import weave

import glob
from datetime import datetime, timedelta
from math import ceil
from itertools import izip

from computeQuantities import theta2Temperature, computeReflectivity, toRecArray

def loadAndComputeCovariance(obs, ens1, ens2, map, axes, normalize):
    cov_shape = [ obs.shape[0] ]
    cov_shape.extend(ens1.shape[1:])

    cov = np.empty(tuple(cov_shape), dtype=np.float32)

    obs_x, obs_y = map(obs['longitude'], obs['latitude'])

    covariance_code = open("cov.c", 'r').read()
    pylib_code = open("pylib.c", 'r').read()

    ens_obs_shape = ( ens1.shape[0], obs.shape[0] )
    ens_obs = np.empty(ens_obs_shape, dtype=np.float32)
    for lde in xrange(ens1.shape[0]):
        ens_obs[lde] = interpolate(ens1[lde], axes, { 'z':obs['elevation'], 'y':obs_y, 'x':obs_x }, wrap=True)

    ens_data = ens2.astype(np.float32)
    for obs_idx in xrange(obs.shape[0]):
        print "Working on ob %d ..." % obs_idx
        ens_ob = ens_obs[:, obs_idx]
        ens_cov = cov[obs_idx]

        weave.inline("ens_covariance(ens_data_array, ens_ob_array, ens_cov_array, PyObject_IsTrue(normalize) == 1);",
            ['ens_data', 'ens_ob', 'ens_cov', 'normalize'],
            support_code=pylib_code + covariance_code,
            libraries=['m'],
        )

#   last_ob_idx = -1
#   for idx in np.ndindex(cov.shape):
#       if idx[0] != last_ob_idx:
#           print "Working on ob %d ..." % idx[0]
#           last_ob_idx = idx[0]

#       ens_idx = [ slice(None) ]
#       ens_idx.extend(idx[1:])
#       ens_var = ens[ens_idx]

#       cov[idx] = np.cov(np.vstack((ens_obs[:, idx[0]], ens_var)))[0, 1]

    return cov

def plotCov(cov, refl, point, map, grid_spacing, file_name, roi=50, normalize=(-1, 1)):
    pylab.figure(figsize=(10, 8))

    pt_z, pt_lat, pt_lon = point
    pt_x, pt_y = map(pt_lon, pt_lat)

    dx, dy = grid_spacing
    nx, ny = cov.shape

    xs, ys = np.meshgrid(dx * np.arange(nx), dy * np.arange(ny))

    lbound, ubound = normalize
    dbound = (ubound - lbound) / 20.
    levels = np.arange(lbound, ubound + dbound, dbound)

    map.contourf(xs, ys, cov, levels=levels)
    pylab.colorbar()
    pylab.contour(xs, ys, refl, levels=[20., 40.], colors='k')

    pylab.plot(pt_x, pt_y, 'ko')
    roi_circle = CirclePolygon((pt_x, pt_y), roi * 1000., resolution=40, ec='k', fc='none')
    axes_box = Rectangle((0, 0), 1, 1, ec='w', fc='w', alpha=0.5, clip_path=roi_circle, transform=(pylab.gca().transAxes + pylab.gca().transData.inverted()))

    path_codes = [ Path.MOVETO ] + (axes_box.get_verts().shape[0] - 1) * [ Path.LINETO ] + [ Path.MOVETO ] + (roi_circle.get_verts().shape[0] - 1) * [ Path.LINETO ]
    path_verts = np.concatenate((axes_box.get_verts(), roi_circle.get_verts()[::-1]))

    mask_path = Path(path_verts, path_codes)
    mask_patch = PathPatch(mask_path, fc='w', ec='w', alpha=0.7)
    pylab.gca().add_patch(mask_patch)

    map.drawcoastlines(linewidth=1.5)
    map.drawcountries(linewidth=1.5)
    map.drawstates(linewidth=1.0)
    if not hasattr(map, 'counties'):
        map.readshapefile('countyp020', 'counties', linewidth=0.5)
    else:
        for county, data in izip(map.counties_info, map.counties):
            if county['STATE'] in ['NE', 'WY', 'CO']:
                pylab.gca().add_patch(Polygon(data, ec='k', fc='none', linewidth=0.5))

    pylab.savefig(file_name)
    pylab.close()
    return

def main():
    times_sec = range(11100, 14700, 300)
    run_base_time = datetime(2009, 6, 5, 18, 0, 0)
    times_dt = [ run_base_time + timedelta(seconds=t) for t in times_sec ]
    radar_elev, radar_lat, radar_lon = 1883, 41.151944, -104.806111
    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs)

    variables = [ 'pt' ]
    interp_height = 25

    if len(variables) == 2:
        var1, var2 = variables
    else:
        var1 = variables[0]
        var2 = variables[0]

    map = Basemap(**proj)
    radar_x, radar_y = map(radar_lon, radar_lat)

    obs_file_names = ['psu_straka_mesonet.pkl', 'ttu_sticknet.pkl', 'asos.pkl', 'soundings_da.pkl']
    all_obs = loadObs(obs_file_names, times_dt, map, sounding_obs=['soundings_da.pkl'])

    print all_obs

    forecast_base = "/caps1/tsupinie/1km-control-20120712/"
    grdbas_file_name = "/caps1/tsupinie/1km-control-20120712/ena001.hdfgrdbas"

    enf_files = glob.glob("%s/enf???.hdf0*" % forecast_base)

    enf_data, ens_members, ens_times = loadAndInterpolateEnsemble(enf_files, variables, toRecArray, grdbas_file_name, 
        points=None, agl=True, wrap=False)

    refl_ens_mean, ens_refl, ens_members, ens_times = loadAndInterpolateEnsemble(enf_files, ['pt', 'p', 'qr', 'qs', 'qh'], computeReflectivity, "/caps1/tsupinie/1km-control-20120712/ena001.hdfgrdbas", 
        {'z':interp_height}, agl=True, wrap=False, aggregator=lambda x: np.mean(x, axis=0))

    grdbas = nio.open_file(grdbas_file_name, mode='r', format='hdf')
    axes_msl = { 'z':grdbas.variables['zp'][:], 'y':grdbas.variables['y'][:], 'x':grdbas.variables['x'][:] }
    axes_agl = { 'z':_makeZCoordsAGL(grdbas.variables['zp'][:]), 'y':grdbas.variables['y'][:], 'x':grdbas.variables['x'][:] }

    for wdt, time_dt in enumerate(times_dt):
        print "Working on time %s" % str(time_dt)
        time_idxs = np.where(all_obs['time'] == (time_dt - datetime(1970, 1, 1, 0, 0, 0)).total_seconds())
        cov = loadAndComputeCovariance(all_obs[time_idxs], enf_data[var1][:, wdt], enf_data[var2][:, wdt], map, axes_msl, False)
        cor = loadAndComputeCovariance(all_obs[time_idxs], enf_data[var1][:, wdt], enf_data[var2][:, wdt], map, axes_msl, True)

        for ob_idx in xrange(cov.shape[0]):
            ob = all_obs[time_idxs[0][ob_idx]]
            interp_cov = interpolate(cov[ob_idx], axes_agl, { 'z':interp_height })
            interp_cor = interpolate(cor[ob_idx], axes_agl, { 'z':interp_height })

            roi = 50
            if ob['obtype'] == "SA":
                roi = 300
            elif ob['obtype'] == "SNDG":
                obs_x, obs_y = map(ob['longitude'], ob['latitude'])
                interp_height_msl = interpolate(axes_msl['z'], axes_msl, { 'y':obs_y, 'x':obs_x }, wrap=True)[0]
                print interp_height_msl

                r_h = 150.
                r_v = 6.

                z_diff = (interp_height_msl - ob['elevation']) / 1000.
                if np.abs(z_diff) > r_v:
                    roi = 0
                else:
                    roi = r_h * np.sqrt(1 - (z_diff / r_v) ** 2)

                print roi
            print np.nanmin(interp_cov), np.nanmax(interp_cov)

            plotCov(interp_cov, refl_ens_mean[wdt], (ob['elevation'], ob['latitude'], ob['longitude']), map, goshen_1km_gs, "cov_%06d_ob%02d.png" % (times_sec[wdt], ob_idx), roi=roi, normalize=(-2.0, 2.0))
            plotCov(interp_cor, refl_ens_mean[wdt], (ob['elevation'], ob['latitude'], ob['longitude']), map, goshen_1km_gs, "cor_%06d_ob%02d.png" % (times_sec[wdt], ob_idx), roi=roi)

    return

if __name__ == "__main__":
    main()
