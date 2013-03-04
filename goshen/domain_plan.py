
import numpy as np

import Nio as nio

import matplotlib
matplotlib.use('agg')
import pylab
from matplotlib.patches import Circle

from mpl_toolkits.basemap import Basemap

import struct

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

def createBasemap(grid_spacing, grid_number, center_lat, center_lon, true_lon):
    domain = Basemap(projection='lcc', resolution=None,
        width=(grid_number * grid_spacing), height=(grid_number * grid_spacing), lat_0=center_lat, lon_0=true_lon,
        lat_1=30., lat_2=60.)

    ctr_x, ctr_y = domain(center_lon, center_lat)
    move_x = ctr_x - grid_number * grid_spacing / 2
    move_y = ctr_y - grid_number * grid_spacing / 2

    llcrnrlon, llcrnrlat = domain(move_x, move_y, inverse=True)
    urcrnrlon, urcrnrlat = domain(grid_number * grid_spacing + move_x, grid_number * grid_spacing + move_y, inverse=True)

    domain = Basemap(projection='lcc', resolution=None,
        llcrnrlat=llcrnrlat, llcrnrlon=llcrnrlon, urcrnrlat=urcrnrlat, urcrnrlon=urcrnrlon,
        lat_0=center_lat, lon_0=true_lon, lat_1=30., lat_2=60.) 

    return domain

def createBoundary(grid_number, lb_x=0, lb_y=0):
    boundary_x = lb_x * np.ones((4 * grid_number,))
    boundary_y = lb_y * np.ones((4 * grid_number,))

    boundary_x[:grid_number] = np.arange(lb_x, lb_x + grid_number)
    boundary_x[grid_number:(2 * grid_number)] = lb_x + grid_number
    boundary_x[(2 * grid_number):(3 * grid_number)] = np.arange(lb_x + grid_number, lb_x, -1)

    boundary_y[grid_number:(2 * grid_number)] = np.arange(lb_y, lb_y + grid_number)
    boundary_y[(2 * grid_number):(3 * grid_number)] = lb_y + grid_number
    boundary_y[(3 * grid_number):(4 * grid_number)] = np.arange(lb_y + grid_number, lb_y, -1)

    return boundary_x, boundary_y

def computeBoundaryCoords(grid_spacing, grid_number, center_lat, center_lon, true_lon, subset=None):
    if subset is not None:
        bounds_x, bounds_y = subset
        subset_grid_number = bounds_x.stop - bounds_x.start   
        boundary_x, boundary_y = createBoundary(subset_grid_number, bounds_x.start, bounds_y.start)
    else:
        boundary_x, boundary_y = createBoundary(grid_number)    

    domain = createBasemap(grid_spacing, grid_number, center_lat, center_lon, true_lon)

    boundary_lons, boundary_lats = domain(grid_spacing * boundary_x, grid_spacing * boundary_y, inverse=True)
    return boundary_lats, boundary_lons

def plotMap(analysis_boundary, inner_boundary, outer_boundary, map_params, radars, radar_range, file_name, grids=None, topo=None):
    pylab.clf()
    pylab.axes((0, 0, 1, 1))
    analysis_boundary_lats, analysis_boundary_lons = analysis_boundary
    inner_boundary_lats, inner_boundary_lons = inner_boundary
    outer_boundary_lats, outer_boundary_lons = outer_boundary
    edge, ctr_lat, ctr_lon = map_params

    map = Basemap(projection='lcc', resolution='i', area_thresh=10000,
        width=edge, height=edge,
        lat_0=ctr_lat, lon_0=ctr_lon, lat_1=30., lat_2=60.)

    analysis_map_x, analysis_map_y = map(analysis_boundary_lons, analysis_boundary_lats)
    inner_map_x, inner_map_y = map(inner_boundary_lons, inner_boundary_lats)
    outer_map_x, outer_map_y = map(outer_boundary_lons, outer_boundary_lats)

    if topo is not None:
        topo_data, topo_lats, topo_lons = topo
        topo_lats, topo_lons = np.meshgrid(topo_lats, topo_lons)
        xs, ys = map(topo_lons, topo_lats)
        pylab.contourf(xs, ys, topo_data, cmap=matplotlib.cm.get_cmap('gist_earth'))
#       pylab.colorbar()

    for radar_id, (radar_lat, radar_lon) in radars.iteritems():
        radar_x, radar_y = map(radar_lon, radar_lat)
        map.plot(radar_x, radar_y, 'ko')
        if radar_id == "MWR-05XP":
            pylab.text(radar_x + 10000, radar_y + 10000, radar_id, ha='left', va='bottom')
        else:
            pylab.text(radar_x - 10000, radar_y - 10000, radar_id, ha='right', va='top')
        pylab.gca().add_patch(Circle((radar_x, radar_y), radius=(radar_range[radar_id] * 1000), edgecolor='k', fill=False)) 

    if grids is not None:
        colors = [ 'g', 'm', 'c'] 
        for (grid_data, grid_lats, grid_lons), color in zip(grids, colors[:len(grids)]):
            grid_x, grid_y = map(grid_lons, grid_lats)
            map.contour(grid_x, grid_y, grid_data, colors=color, levels=range(30, 80, 10))

    map.plot(analysis_map_x, analysis_map_y, "#333333")
    map.plot(inner_map_x, inner_map_y, 'r')
    map.plot(outer_map_x, outer_map_y, 'b')

    map.drawcoastlines()
    map.drawcountries()
    map.drawstates()
#   map.drawparallels(np.arange(20, 60, 5), dashes=[2,1])
#   map.drawmeridians(np.arange(-120, -90, 5), dashes=[2,1])

    pylab.savefig(file_name)
    return

def computeGridCoords(grid_spacing, grid_number, center_lat, center_lon, offset_x=0, offset_y=0):
    grid_number_x, grid_number_y = grid_number
    domain = Basemap(projection='lcc', resolution=None,
        width=(grid_number_x * grid_spacing), height=(grid_number_y * grid_spacing), lat_0=center_lat, lon_0=center_lon,
        lat_1=30., lat_2=60.)

    coords_x, coords_y = np.meshgrid(np.arange(grid_number_x) - offset_x, np.arange(grid_number_y) - offset_y)

#   print grid_number
#   print np.arange(grid_number_x) - offset_x
#   print np.arange(grid_number_y) - offset_y

    lons, lats = domain(grid_spacing * coords_x, grid_spacing * coords_y, inverse=True)
    return lats, lons

def loadGrid(file_name, center_lat, center_lon, subset=None):
    hdf = nio.open_file(file_name, mode='r', format='hdf')
    grid = hdf.variables['refl2d'][0]

    if subset is None:
        offset_x = 0
        offset_y = 0
    else:
        full_grid_y, full_grid_x = grid.shape

        grid = grid[subset[1], subset[0]]

        offset_x = full_grid_x / 2 - (subset[0].start + subset[0].stop) / 2
        offset_y = full_grid_y / 2 - (subset[1].start + subset[1].stop) / 2

    lats, lons = computeGridCoords(1000., tuple(reversed(grid.shape)), center_lat, center_lon, offset_x=offset_x, offset_y=offset_y)

    return grid, lats, lons

def main():
    
    radar_location = {'KCYS':(41.15194, -104.80611), 'KFTG':(39.78667, -104.54583), 'KRIW':(43.06611, -108.47722), 'MWR-05XP':(41.56150, -104.298996)} #, 'KUDX':(44.125, -102.82972)}
    radar_range = {'KCYS':230, 'KFTG':230, 'KRIW':230, 'MWR-05XP':40}

    topo, topo_lats, topo_lons = load_topo("e10g", (6000, 10800), ((0., 50.), (-180., -90.)), ((30., 50.), (-120., -94.)))

    analysis_subset = (slice(90, 170), slice(100, 180))

    inner_spacing = 1000
    number_inner = 252
    inner_center_lat = 41.61795
    inner_center_lon = -104.34843

    outer_spacing = 3000
    number_outer = 408
    outer_center_lat = 40.61999
    outer_center_lon = -107.344

    subset_ll_x = 110
    subset_ll_y = 150
    subset_ur_x = 200
    subset_ur_y = 250

#   initial_grid = loadGrid("goshen.hdfrefl2d000000", inner_center_lat, inner_center_lon, subset=(slice(subset_ll_x, subset_ur_x), slice(subset_ll_y, subset_ur_y)))

    subset_ll_x = 200
    subset_ll_y = 175
    subset_ur_x = 316
    subset_ur_y = 227

#   final_grid = loadGrid("goshen.hdfrefl2d014818", inner_center_lat, inner_center_lon,  subset=(slice(subset_ll_x, subset_ur_x), slice(subset_ll_y, subset_ur_y)))

    analysis_boundary = computeBoundaryCoords(inner_spacing, number_inner, inner_center_lat, inner_center_lon, outer_center_lon, subset=analysis_subset)
    inner_boundary = computeBoundaryCoords(inner_spacing, number_inner, inner_center_lat, inner_center_lon, outer_center_lon)
    outer_boundary = computeBoundaryCoords(outer_spacing, number_outer, outer_center_lat, outer_center_lon, outer_center_lon)

    plotMap(analysis_boundary, inner_boundary, outer_boundary, (outer_spacing * (100 + number_outer), outer_center_lat, outer_center_lon), radar_location, radar_range, "domain_plan.png", topo=(topo, topo_lats, topo_lons))
    plotMap(analysis_boundary, inner_boundary, outer_boundary, (inner_spacing * (100 + number_inner), inner_center_lat, inner_center_lon), radar_location, radar_range, "domain_plan_zoom.png") #,
         #grids=[initial_grid, final_grid])

    return

if __name__ == "__main__":
    main()
