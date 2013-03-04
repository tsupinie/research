
import struct

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

from mpl_toolkits.basemap import Basemap

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

    keep_2d_jdys, keep_2d_idxs = np.meshgrid(keep_jdys[0], keep_idxs[0])

    return topo_data[keep_2d_jdys, keep_2d_idxs], lats[keep_jdys], lons[keep_idxs]

if __name__ == "__main__":
    topo_data, lats, lons = load_topo("e10g", (6000, 10800), ((0., 50.), (-180., -90.)), ((36., 46.), (-109., -99.)))
    lats, lons = np.meshgrid(lats, lons)

    print topo_data.shape
    print lats.shape
    print lons.shape

    map = Basemap(projection='merc', resolution='l', llcrnrlat=lats.min(), llcrnrlon=lons.min(), urcrnrlat=lats.max(), urcrnrlon=lons.max(), lat_ts=0.)

    topo_x, topo_y = map(lons, lats)
    pylab.contourf(topo_x, topo_y, topo_data)

    map.drawcoastlines()
    map.drawcountries()
    map.drawstates()
    pylab.savefig("topo_test.png")
