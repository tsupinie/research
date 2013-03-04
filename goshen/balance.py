
import numpy as np

import Nio as nio

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap

from glob import glob
import struct

def verticalDerivative(data, height):
    derivative = np.zeros(data.shape)
    alpha = (height[2:] - height[1:-1]) / (height[1:-1] - height[:-2])

    derivative[1:-1] = (data[2:] + (alpha - 1) * data[1:-1] - alpha * data[:-2]) / (2 * alpha * (height[1:-1] - height[:-2]))
    derivative[0] = (data[1] - data[0]) / (height[1] - height[0])
    derivative[-1] = (data[-1] - data[-2]) / (height[-1] - height[-2])

    return derivative

def computeBalances(data_file, grdbas_file):
    # Balances: thermal wind, gradient wind, hydrostatic, mass conservation

    gas_const_dry = 287.
    accel_gravity = 9.806
    coriolis_param = 2 * (2 * np.pi / 86400) * np.sin(np.pi * 40. / 180.)

    temperature = data_file.variables['pt'][:,:-1,:-1] * (data_file.variables['p'][:,:-1,:-1] / 100000.) ** (2. / 7.)
    density = data_file.variables['p'][:,:-1,:-1] / (gas_const_dry * temperature)

    hydrostatic_lhs = verticalDerivative(data_file.variables['p'][:,:-1,:-1], grdbas_file.variables['zp'][:,:-1,:-1])
    hydrostatic_rhs = -density * accel_gravity

    hydrostatic = hydrostatic_lhs - hydrostatic_rhs    

    dx = grdbas_file.variables['x'][1] - grdbas_file.variables['x'][0]
    dy = grdbas_file.variables['y'][1] - grdbas_file.variables['y'][0]

    # Incompressible mass conservation balance
    mass_cons_lhs = verticalDerivative(data_file.variables['w'][:,:-1,:-1], grdbas_file.variables['zp'][:,:-1,:-1])
    mass_cons_rhs = -np.gradient(data_file.variables['u'][:,:-1,:-1], dx)[1] - np.gradient(data_file.variables['v'][:,:-1,:-1], dy)[2]

    mass_cons = mass_cons_lhs - mass_cons_rhs

    # Thermal wind balance
    dummy, dp_dx, dp_dy = np.gradient(data_file.variables['p'][:,:-1,:-1], 1, dx, dy)
    dummy, dT_dx, dT_dy = np.gradient(temperature, 1, dx, dy)
    u_g = -1 / (coriolis_param * density) * dp_dy
    v_g = 1 / (coriolis_param * density) * dp_dx

    thermal_wind_u_lhs = verticalDerivative(u_g, grdbas_file.variables['zp'][:,:-1,:-1])
    thermal_wind_u_rhs = -(accel_gravity * dT_dy) / (coriolis_param * temperature) + u_g / temperature * verticalDerivative(temperature, grdbas_file.variables['zp'][:,:-1,:-1])

    thermal_wind_u = thermal_wind_u_lhs - thermal_wind_u_rhs

    thermal_wind_v_lhs = verticalDerivative(v_g, grdbas_file.variables['zp'][:,:-1,:-1])
    thermal_wind_v_rhs = (accel_gravity * dT_dx) / (coriolis_param * temperature) + v_g / temperature * verticalDerivative(temperature, grdbas_file.variables['zp'][:,:-1,:-1])

    thermal_wind_v = thermal_wind_v_lhs - thermal_wind_v_rhs

    return hydrostatic, mass_cons, thermal_wind_u, thermal_wind_v

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

def plot_map(data_grid, grid_spacing, orientation, title, file_name, topo=None):
    pylab.clf()
    nx, ny = data_grid.shape
    grid_x, grid_y = grid_spacing

    map = Basemap(projection='lcc', resolution=None, width=(nx * grid_x), height=(ny * grid_y),
        lat_0=40.61998, lon_0=-107.344, lat_1=30., lat_2=60.) 

    ctr_x, ctr_y = map(-104.344, 41.61998)
    move_x = ctr_x - nx * grid_x / 2
    move_y = ctr_y - ny * grid_y / 2

    llcrnrlon, llcrnrlat = map(move_x, move_y, inverse=True)
    urcrnrlon, urcrnrlat = map(nx * grid_x + move_x, ny * grid_y + move_y, inverse=True)

    map = Basemap(projection='lcc', resolution='l',
        llcrnrlat=llcrnrlat, llcrnrlon=llcrnrlon, urcrnrlat=urcrnrlat, urcrnrlon=urcrnrlon,
        lat_0=40.61998, lon_0=-107.344, lat_1=30., lat_2=60.) 

    x, y = np.meshgrid(grid_x * np.arange(nx), grid_y * np.arange(ny))
    interval = (data_grid.max() - data_grid.min()) / 100.
    pylab.contourf(x, y, data_grid, levels=np.arange(data_grid.min(), data_grid.max() + interval, interval))

    if topo is not None:
        topo_data, topo_lats, topo_lons = topo

        topo_lats, topo_lons = np.meshgrid(topo_lats, topo_lons)
        topo_x, topo_y = map(topo_lons, topo_lats)

        map.contour(topo_x, topo_y, topo_data, colors='#808080')

    map.drawstates()
    map.drawcountries()
    map.drawcoastlines()

    pylab.colorbar()
    pylab.title(title)
    pylab.savefig(file_name)
    return

def main():
    files = glob("1kmgoshen/1kmgoshen.hdf0*")
    hdf_grdbas = nio.open_file("1kmgoshen/1kmgoshen.hdfgrdbas", mode='r', format='hdf')

#   topo, topo_lats, topo_lons = load_topo("e10g", (6000, 10800), ((0., 50.), (-180., -90.)), ((36., 46.), (-114., -101.)))

    for file in files:
        time_sec = file[-6:]
        hdf_data = nio.open_file(file, mode='r', format='hdf')
        hydrostatic, mass_cons, thermal_wind_u, thermal_wind_v = computeBalances(hdf_data, hdf_grdbas)
        plot_map(hydrostatic[1], (1000, 1000), "xy", r"Hydrostatic Imbalance (Pa m$^{-1}$)", "hydrostatic_t%s.png" % time_sec) #, topo=(topo, topo_lats, topo_lons))
        plot_map(mass_cons[1], (1000, 1000), "xy", r"Mass Conservation Imbalance (m s$^{-2}$)", "mass_cons_t%s.png" % time_sec) #, topo=(topo, topo_lats, topo_lons))
        plot_map(thermal_wind_u[1], (1000, 1000), "xy", r"Thermal Wind $u$ Imbalance (m s$^{-2}$)", "thermal_wind_u_t%s.png" % time_sec) #, topo=(topo, topo_lats, topo_lons))
        plot_map(thermal_wind_v[1], (1000, 1000), "xy", r"Thermal Wind $v$ Imbalance (m s$^{-2}$)", "thermal_wind_v_t%s.png" % time_sec) #, topo=(topo, topo_lats, topo_lons))
        hdf_data.close()

    hdf_grdbas.close()
    return

if __name__ == "__main__":
    main()
