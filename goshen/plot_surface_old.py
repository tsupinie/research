import Nio as nio

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap

files = {
    "3km-ensemble-20120511/enmean.hdf010800": 'enmean-20120511.sfc.010800.png', 
    "3km-ensemble-20120531/enmean.hdf010800": 'enmean-20120531.sfc.010800.png', 
    "3km-ensemble-20120608/enmean.hdf010800": 'enmean-20120608.sfc.010800.png',
    "3kmgoshenens.hdf000000": 'arps-anal.sfc.000000.png',
    "3kmgoshenicbc.hdf021600": 'nam-anal.sfc.000000.png',
}

for input_file, output_file in files.iteritems():
    hdf = nio.open_file(input_file, mode='r', format='hdf')

    u = hdf.variables['u'][1]
    v = hdf.variables['v'][1]
    p = hdf.variables['p'][1]
    qv = hdf.variables['qv'][1]

#   print "u", np.unravel_index(np.argmax(u), u.shape), u.max()
#   print "v", np.unravel_index(np.argmax(v), v.shape), v.max()

#   u = u[1]
#   v = v[1]

#   V = np.sqrt(u ** 2 + v ** 2)

    stride = 15
    grid_spacing = 3000.
    nx, ny = u.shape

    x = np.arange(nx) * grid_spacing
    y = np.arange(ny) * grid_spacing

    x, y = np.meshgrid(x, y)

    map = Basemap(projection='lcc', resolution='l', width=(nx * grid_spacing), height=(ny * grid_spacing),
        lat_0=40.61998, lon_0=-107.344, lat_1=30., lat_2=60.) 

#   ctr_x, ctr_y = map(-104.34843, 41.61975)
#   move_x = ctr_x - nx * grid_spacing / 2
#   move_y = ctr_y - ny * grid_spacing / 2

#   llcrnrlon, llcrnrlat = map(move_x, move_y, inverse=True)
#   urcrnrlon, urcrnrlat = map(nx * grid_spacing + move_x, ny * grid_spacing + move_y, inverse=True)

#   map = Basemap(projection='lcc', resolution='l',
#       llcrnrlat=llcrnrlat, llcrnrlon=llcrnrlon, urcrnrlat=urcrnrlat, urcrnrlon=urcrnrlon,
#       lat_0=40.61998, lon_0=-107.344, lat_1=30., lat_2=60.) 

    pylab.clf()
    map.contourf(x, y, qv, levels=np.arange(0, 0.0135, 0.0015))
    pylab.colorbar()
    map.contour(x, y, p, colors='k')
    map.quiver(x[::stride, ::stride], y[::stride, ::stride], u[::stride, ::stride], v[::stride, ::stride])
    map.drawstates()
    pylab.savefig(output_file)

    hdf.close()
