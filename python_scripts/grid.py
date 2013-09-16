
def _imports():
    global copy, izip, np, pylab, Polygon, Basemap

    import numpy as np
    import pylab
    from matplotlib.patches import Polygon
    from mpl_toolkits.basemap import Basemap
    
    from itertools import izip
    import copy
    return

if __name__ != "__main__":
    _imports()

class ExperimentGrid:
    def __init__(self, **kwargs):
        if 'bounds' in kwargs and kwargs['bounds']:
            self._bound_x, self._bound_y = kwargs['bounds']
        else:
            self._bound_x, self._bound_y = slice(None), slice(None)

        if 'bounds' in kwargs:
            del kwargs['bounds']

        self._nx, self._ny = kwargs['grid_size']
        self._dx, self._dy = kwargs['grid_spacing']
        del kwargs['grid_size']
        del kwargs['grid_spacing']

        kwargs['rsphere'] = 6371000

        kwargs['width'] = self._dx * self._nx
        kwargs['height'] = self._dy * self._ny

        if 'ctr_lat' in kwargs:
            self._ctr_lat = kwargs['ctr_lat']
            del kwargs['ctr_lat']
        else:
            self._ctr_lat = kwargs['lat_0']

        if 'ctr_lon' in kwargs:
            self._ctr_lon = kwargs['ctr_lon']
            del kwargs['ctr_lon']
        else:
            self._ctr_lon = kwargs['lon_0']

        if 'state_list' in kwargs:
            self._state_list = kwargs['state_list']
            del kwargs['state_list']
        else:
            self._state_list = None

        self._proj = copy.copy(kwargs)
        map = Basemap(**kwargs) 

        ctr_x, ctr_y = map(self._ctr_lon, self._ctr_lat)
        move_x = ctr_x - kwargs['width'] / 2
        move_y = ctr_y - kwargs['height'] / 2

        def defaultValue(value, default):
            if value is not None:  return value
            else:                  return default

        lower_bound_x = defaultValue(self._bound_x.start, 0)
        lower_bound_y = defaultValue(self._bound_y.start, 0)
        upper_bound_x = defaultValue(self._bound_x.stop, self._nx) - 1
        upper_bound_y = defaultValue(self._bound_y.stop, self._ny) - 1

        llcrnrlon, llcrnrlat = map(move_x + self._dx * lower_bound_x, move_y + self._dy * lower_bound_y, inverse=True)
        urcrnrlon, urcrnrlat = map(move_x + self._dx * upper_bound_x, move_y + self._dy * upper_bound_y, inverse=True)

        del kwargs['width']
        del kwargs['height']

        kwargs['llcrnrlat'] = llcrnrlat
        kwargs['llcrnrlon'] = llcrnrlon
        kwargs['urcrnrlat'] = urcrnrlat
        kwargs['urcrnrlon'] = urcrnrlon

        self._map = Basemap(**kwargs)
        return

    def getXY(self, x=None, y=None):
        xs, ys = np.meshgrid(self._dx * (np.arange(self._nx) - 1), self._dy * (np.arange(self._ny) - 1))

        xs = xs[self.getBounds()]
        xs = xs - xs[0, 0]

        ys = ys[self.getBounds()]
        ys = ys - ys[0, 0]

        if x is None or y is None:
            return xs, ys
        else:
            return xs[x, y], ys[x, y]

    def getBounds(self):
        return self._bound_y, self._bound_x

    def getBoundaryCoords(self):
        def _createBoundary(nx, ny, lb_x=0, lb_y=0):
            boundary_x = lb_x * np.ones((2 * (nx + ny),))
            boundary_y = lb_y * np.ones((2 * (nx + ny),))

            boundary_x[:nx] = np.arange(lb_x, lb_x + nx)
            boundary_x[nx:(nx + ny)] = lb_x + nx 
            boundary_x[(nx + ny):(2 * nx + ny)] = np.arange(lb_x + nx, lb_x, -1)

            boundary_y[nx:(nx + ny)] = np.arange(lb_y, lb_y + ny)
            boundary_y[(nx + ny):(2 * nx + ny)] = lb_y + ny
            boundary_y[(2 * nx + ny):(2 * (nx + ny))] = np.arange(lb_y + ny, lb_y, -1)

            return boundary_x, boundary_y

        if self._bound_x != slice(None) and self._bound_y != slice(None):
            subset_nx = self._bound_x.stop - self._bound_x.start   
            subset_ny = self._bound_y.stop - self._bound_y.start   
            boundary_x, boundary_y = _createBoundary(subset_nx, subset_ny, self._bounds_x.start, self._bound_y.start)
        else:
            boundary_x, boundary_y = _createBoundary(self._nx, self._ny)

        boundary_lons, boundary_lats = self(self._dx * boundary_x, self._dy * boundary_y, inverse=True)
        return boundary_lats, boundary_lons

    def getGridSpacing(self):
        return self._dx, self._dy

    def getWidthHeight(self, override=False):
        if self._bound_x == slice(None) and self._bound_y == slice(None) or override:
            return self._dx * self._nx, self._dy * self._ny
        else:
            return self._dx * (self._bound_x.stop - self._bound_x.start), self._dy * (self._bound_y.stop - self._bound_y.start)

    def drawPolitical(self, scale_len=0, color='k', lw=1):
        self._map.drawcountries(linewidth=1.0, color=color)
        self._map.drawcoastlines(linewidth=1.0, color=color)

        if not hasattr(self._map, 'counties'):
            self._map.readshapefile("countyp020", 'counties', linewidth=(lw * 0.5), color=color)
        else:
            for county, data in izip(self._map.counties_info, self._map.counties):
                if not self._state_list or county['STATE'] in self._state_list:
                    pylab.gca().add_patch(Polygon(data, ec=color, fc='none', linewidth=(lw * 0.5)))

        if not hasattr(self._map, 'states'):
            self._map.readshapefile("tl_2011_us_state", 'states', linewidth=(lw * 1.5), color=color)
        else:
            for state, data in izip(self._map.states_info, self._map.states):
                if not self._state_list or state['STUSPS'] in self._state_list:
                    pylab.gca().add_patch(Polygon(data, ec=color, fc='none', linewidth=(lw * 1.5)))

        if scale_len > 0:
            offset_x, offset_y = 0.03, 0.09
            ax_trans = pylab.gca().transData + pylab.gca().transAxes.inverted() 
            coords_x, coords_y = ax_trans.transform(1000 * np.array([[0, scale_len], [0, 0]]))
            pylab.plot(coords_x + offset_x, coords_y + offset_y, 'k-', lw=(lw * 2.5), transform=pylab.gca().transAxes)
            pylab.text(coords_x.mean() + offset_x, coords_y.mean() + offset_y - 0.025, "%d km" % scale_len, ha='center', va='top', transform=pylab.gca().transAxes)
        return

    def __call__(self, *args, **kwargs):
        if 'inverse' in kwargs and kwargs['inverse']:
            return self._map(*[ a + d for a, d in zip(args, [1.5 * self._dx, 1.5 * self._dy]) ], **kwargs)
        else:
            return [ c - d for c, d in zip(self._map(*args, **kwargs), [1.5 * self._dx, 1.5 * self._dy]) ]

goshen_1km_grid = lambda bounds=None, buffer=0: \
    ExperimentGrid(projection='lcc', resolution='l',
    grid_size=(255 + 2 * buffer, 255 + 2 * buffer), grid_spacing=(1000, 1000),
    ctr_lat=41.61975, ctr_lon=-104.34843, 
    lat_0=40.619985, lon_0=-107.344, lat_1=30., lat_2=60.,
    state_list=['CO', 'NE', 'WY'], bounds=bounds)

goshen_3km_grid = lambda bounds=None, buffer=0: \
    ExperimentGrid(projection='lcc', resolution='l',
    grid_size=(411 + 2 * buffer, 411 + 2 * buffer), grid_spacing=(3000, 3000),
    lat_0=40.619985, lon_0=-107.344, lat_1=30., lat_2=60.,
    state_list=['CO', 'NE', 'WY', 'SD', 'MT', 'ID', 'UT', 'NV', 'AZ', 'NM', 'TX', 'OK', 'KS' ], bounds=bounds)

if __name__ == "__main__":
    import matplotlib
    matplotlib.use('agg')
    _imports()

    grid = goshen_1km_grid()
    grid.drawPolitical()
    pylab.savefig("grid_test.png")
