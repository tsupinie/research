
import glob
from datetime import datetime, timedelta
from multiprocessing import Process
import struct
import csv
import os.path

import Nio as nio

import numpy as np

from computeQuantities import computeReflectivity
from legacy import drawPolitical, goshen_1km_proj, goshen_1km_gs, goshen_3km_proj, goshen_3km_gs, setupMapProjection
from util import decompressVariable
from color_tables import NWSRef

import pylab
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Circle

from StationReader import StationReader

def load_topo(file_name, grid_size, bounds, trim_bounds):
    nx, ny = grid_size
    lat_bounds, lon_bounds = bounds
    lat_trim, lon_trim = trim_bounds

    file = open(file_name, 'r')
    topo_string = file.read()
    file.close()

    topo_data = np.array(struct.unpack('<' + 'h' * (len(topo_string) / 2), topo_string)).reshape((nx, ny))

    lat_lbound, lat_ubound = lat_bounds
    lats = lat_lbound + (np.arange(nx, 0, -1, dtype=np.float32) - 1) / nx * (lat_ubound - lat_lbound)

    lon_lbound, lon_ubound = lon_bounds
    lons = lon_lbound + np.arange(ny, dtype=np.float32) / ny * (lon_ubound - lon_lbound)

    lat_lbound, lat_ubound = lat_trim
    keep_jdys = np.where((lats >= lat_lbound) & (lats <= lat_ubound))

    lon_lbound, lon_ubound = lon_trim
    keep_idxs = np.where((lons >= lon_lbound) & (lons <= lon_ubound))

    return topo_data[np.meshgrid(keep_jdys[0], keep_idxs[0])], lats[keep_jdys], lons[keep_idxs]

def plot_map(radar_data, grid_spacing, title, file_name, color_bar='refl', topo=None, aux_field=None, vectors=None, obs=None):
#   bounds = (slice(80, 160), slice(90, 170))
    bounds = (slice(None), slice(None))
#   bounds = (slice(242, 327), slice(199, 284))
    pylab.figure()
    pylab.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.02)
    nx, ny = radar_data[bounds[::-1]].shape

    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs, bounds=bounds)
    map = Basemap(**proj)

    radar_x, radar_y = map([-104.299004], [41.561497])

    if topo is not None:
        topo_data, topo_lats, topo_lons = topo

        topo_lats, topo_lons = np.meshgrid(topo_lats, topo_lons)
        topo_x, topo_y = map(topo_lons, topo_lats)

        map.contourf(topo_x, topo_y, topo_data, cmap=pylab.get_cmap('gray'))

    if color_bar == 'refl':
        levels = range(10, 85, 5)
        color_map = NWSRef #pylab.get_cmap('jet')
    elif color_bar == 'radv':
        levels = range(-35, 40, 5)
        color_map = pylab.get_cmap('RdBu')
    elif color_bar == 'pt':
        levels = range(296, 321, 2)
        color_map = pylab.get_cmap('jet')

    print nx, ny

    x, y = np.meshgrid(grid_spacing * np.arange(nx), grid_spacing * np.arange(ny))

    map.contourf(x, y, radar_data[bounds[::-1]], levels=levels, cmap=color_map)
#   map.plot(radar_x, radar_y, 'ko')
    pylab.colorbar()

    if aux_field is not None:
        aux_levels, aux_data = aux_field
        CS = map.contour(x, y, aux_data[bounds[::-1]], levels=aux_levels, colors='k', lw=0.75)
#       pylab.clabel(CS, fmt='%.0f', inline_spacing=1, fontsize="12px")

    drawPolitical(map, scale_len=25)

    if obs is not None:
        for stid, ob in obs.iteritems():
            potential_temperature = (5. / 9. * (ob['temperature'] - 32) + 273.15) * (29.5250192 / ob['pressure']) ** (2. / 7.)
            ob_x, ob_y = map(ob['Longitude'], ob['Latitude'])

            ob_ax_x, ob_ax_y = (pylab.gca().transData + pylab.gca().transAxes.inverted()).transform(np.array([ob_x, ob_y]))

            if ob_ax_x > 0 and ob_ax_x <= 1 and ob_ax_y > 0 and ob_ax_y <= 1:
                pylab.gca().add_patch(Circle((ob_x, ob_y), 4000, fc='none', ec='k'))
                pylab.text(ob_x, ob_y, "%5.1f" % potential_temperature, size='x-small', ha='right', va='bottom')

    if vectors is not None:
        stride = 4
        u, v = vectors
        pylab.quiver(x[::stride, ::stride], y[::stride, ::stride], u[::stride, ::stride], v[::stride, ::stride])

    pylab.title(title)
    pylab.savefig(file_name)
    return

def do_observed_plot(file_name, radar_id="KCYS", param="refl", volumetric=False):
    sec_string = file_name[-6:]

    init_time = datetime(2009, 6, 5, 18, 0, 0)
    valid_time = init_time + timedelta(seconds=int(sec_string))

    hdf = nio.open_file(file_name, mode='r', format='hdf')
    if volumetric:
        radar_data = hdf.variables["%s3d" % param][:]
        radar_data_masked = np.ma.array(radar_data, mask=(radar_data == -999))
        hdf.close()

        for idx in range(10, 51):
            if param == 'refl':
                title = "%s Level %d Reflectivity Valid %s" % (radar_id, idx + 1, valid_time.strftime("%d %b %Y %H%M UTC"))
                img_file_name = "%s.l%02dref.%s.png" % (radar_id, idx + 1, sec_string)
            elif param == 'radv':
                title = "%s Level %d Velocity Valid %s" % (radar_id, idx + 1, valid_time.strftime("%d %b %Y %H%M UTC"))
                img_file_name = "%s.l%02dvel.%s.png" % (radar_id, idx + 1, sec_string)

            plot_map(radar_data_masked[idx], 3000, title, img_file_name, color_bar=param)

    else:
        radar_data = hdf.variables["%s2d" % param][0]
        radar_data_masked = np.ma.array(radar_data, mask=(radar_data == -999))
        hdf.close()

        if param == 'refl':
            title = "%s Base Reflectivity Valid %s" % (radar_id, valid_time.strftime("%d %b %Y %H%M UTC"))
            img_file_name = "%s.bref.%s.png" % (radar_id, sec_string)
        elif param == 'radv':
            title = "%s Base Velocity Valid %s" % (radar_id, valid_time.strftime("%d %b %Y %H%M UTC"))
            img_file_name = "%s.bvel.%s.png" % (radar_id, sec_string)

        plot_map(radar_data_masked, 1000, title, img_file_name, color_bar=param)
    return

def load_reflectivity_vars(file_name, vars, ens_member):
    hdf = nio.open_file(file_name, mode='r', format='hdf')

    vars['pt'][ens_member] = hdf.variables['pt'][12]
    vars['p'][ens_member] = hdf.variables['p'][12]
    vars['qr'][ens_member] = np.maximum(hdf.variables['qr'][12], np.zeros(hdf.variables['qr'][12].shape))
    vars['qs'][ens_member] = np.maximum(hdf.variables['qs'][12], np.zeros(hdf.variables['qs'][12].shape))
    vars['qh'][ens_member] = np.maximum(hdf.variables['qh'][12], np.zeros(hdf.variables['qh'][12].shape))

    hdf.close()
    return

def do_simulated_plot(file_name, valid_time, sec_string, obs=None, date_tag=True):
    path_name, file_base = os.path.split(file_name)
    ena_string = file_base.split(".")[0]

    if date_tag:
        tag = path_name.split("/")[-1].split("-")[-1]
        id_string = "-%s" % tag
    else:
        id_string = ""

    hdf = nio.open_file(file_name, mode='r', format='hdf')

    vars = {}
    for var in ['pt', 'p', 'qr', 'qs', 'qh', 'w']:
        vars[var] = hdf.variables[var][12]

        if vars[var].min() == -32768 or vars[var].max() == 32767:
            dindex = (12, slice(None), slice(None))
            vars[var] = decompressVariable(hdf.variables[var], dindex=dindex)

        if var in ['qr', 'qs', 'qh']:
            vars[var] = np.maximum(vars[var], np.zeros(vars[var].shape))

    reflectivity = computeReflectivity(**vars)

    reflectivity_thresh = np.maximum(reflectivity, np.zeros(reflectivity.shape))

    w_levels = range(-20, 0, 2)
    w_levels.extend(range(2, 22, 2))

    refl_title = "Base Reflectivity Valid %s" % valid_time.strftime("%d %b %Y %H%M UTC")
    refl_img_file_name = "%s%s.bref.%s.png" % (ena_string, id_string, sec_string)
    ptprt_title = r"$\theta$ Valid %s" % valid_time.strftime("%d %b %Y %H%M UTC")
    ptprt_img_file_name = "%s%s.pt.%s.png" % (ena_string, id_string, sec_string)

#   vars = {}
    for var in ['pt', 'u', 'v']:
        vars[var] = hdf.variables[var][2]

        if vars[var].min() == -32768 or vars[var].max() == 32767:
            dindex = (2, slice(None), slice(None))
            vars[var] = decompressVariable(hdf.variables[var], dindex=dindex)

#   print pt.min(), pt.max()

    bounds = (slice(100, 130), slice(100, 130))

    plot_map(vars['pt'], 1000, ptprt_title, ptprt_img_file_name, color_bar='pt', vectors=(vars['u'], vars['v']), obs=obs)
    plot_map(reflectivity, 1000, refl_title, refl_img_file_name, color_bar='refl', aux_field=(w_levels, vars['w']))
    return

def main():
    radar_id = "KCYS"
    parameter = "refl"
    volumetric = False
    ensemble_plot = False
    plot_obs = False
#   files = sorted(glob.glob("hdf/%s/1km/goshen.hdf%s2d*" % (radar_id, parameter)))
    files = sorted(glob.glob("/caps2/tsupinie/1kmf-zupdtpt/ena011.hdf*"))
#   files = [ "/data6/tsupinie/adapt=0.80,noise=0.75/enmean.hdf003600" ]#, "/data6/tsupinie/relax=0.50,noise=0.75/enmean.hdf003600" ]

#   topo, topo_lats, topo_lons = load_topo("e10g", (6000, 10800), ((0., 50.), (-180., -90.)), ((36., 46.), (-109., -99.)))

    if plot_obs:
        sfc_stations = StationReader("sfc_stations.csv")
        csv_reader = csv.DictReader(open("obs/1min_obs.csv", 'r'))

        observations = {}
        for row in csv_reader:
            if row['stid'] not in observations:
                observations[row['stid']] = {}

            for key in ['day', 'utc_hour', 'utc_minute', 'temperature', 'dewpoint', '2min_wdir', '2min_wspd']:
                try:
                    observations[row['stid']][key].append(float(row[key]))
                except KeyError:
                    observations[row['stid']][key] = [ float(row[key]) ]

            pressures = [ float(row[b]) for b in ['barometer1', 'barometer2', 'barometer3'] if row[b] != '' ]
            try:
                observations[row['stid']]['pressure'].append(sum(pressures) / len(pressures))
            except KeyError:
                observations[row['stid']]['pressure'] = [ sum(pressures) / len(pressures) ]

            station = sfc_stations.searchByID(row['stid'])

            for key in ['Latitude', 'Longitude', 'Elevation']:
                try:
                    observations[row['stid']][key].append(station[key])
                except:
                    observations[row['stid']][key] = [ station[key] ]

        for stid in observations.keys():
            for param in observations[stid].keys():
                observations[stid][param] = np.array(observations[stid][param])

    procs = []
    max_procs = 8

    if ensemble_plot:
        hdf = nio.open_file(files[0], mode='r', format='hdf')
        grid_shape = hdf.variables['p'].shape
        hdf.close()

        shape = [ len(files) ]
        shape.extend(grid_shape[1:])
        shape = tuple(shape)

        vars = { 'pt':np.empty(shape), 'p':np.empty(shape), 'qr':np.empty(shape), 'qs':np.empty(shape), 'qh':np.empty(shape) }

    for idx, file_name in enumerate(files):
        print "Plotting data in %s ... " % file_name

        sec_string = file_name[-6:]

        init_time = datetime(2009, 6, 5, 18, 0, 0)
        valid_time = init_time + timedelta(seconds=int(sec_string))

        if plot_obs:
            target_day = int(valid_time.strftime("%d"))
            target_hour = int(valid_time.strftime("%H"))
            target_min = int(valid_time.strftime("%M"))

            keep_obs = {}
            for stid, obs in observations.iteritems():
                index = np.where((obs['day'] == target_day) & (obs['utc_hour'] == target_hour) & (obs['utc_minute'] == target_min))[0]
                if index.shape[0] > 0:
                    index = index[0]
                else:
                    continue

                keep_obs[stid] = {}
                for key in obs.keys():
                    keep_obs[stid][key] = obs[key][index]

#       proc = Process(target=do_observed_plot, args=(file_name, radar_id, parameter, volumetric))
        proc = Process(target=do_simulated_plot, args=(file_name, valid_time, sec_string), kwargs={'date_tag':True})#, kwargs={'obs':keep_obs})
        proc.start()
        procs.append(proc)

        if ensemble_plot:
            print "Loading data for mean reflectivity ..."
            load_reflectivity_vars(file_name, vars, idx)

        if len(procs) == max_procs or file_name == files[-1]:
            for proc in procs:
                proc.join()
            procs = []

    if ensemble_plot:
        for key, var in vars.iteritems():
            vars[key] = var.mean(axis=0)

        reflectivity = computeReflectivity(**vars)
        plot_map(reflectivity, 1000, "Ensemble Mean Base Reflectivity Valid 05 June 2009 2300 UTC", '1kmenmean.bref.018000.png')

    return

if __name__ == "__main__":
    main()
