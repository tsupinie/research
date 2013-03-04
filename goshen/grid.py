
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
        upper_bound_x = defaultValue(self._bound_x.stop, self._nx)
        upper_bound_y = defaultValue(self._bound_y.stop, self._ny)

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

    def getXY(self):
        xs, ys = np.meshgrid(self._dx * (np.arange(self._nx) - 1), self._dy * (np.arange(self._ny) - 1))

        xs = xs[self.getBounds()]
        xs = xs - xs[0, 0]

        ys = ys[self.getBounds()]
        ys = ys - ys[0, 0]
        return xs, ys

    def getBounds(self):
        return self._bound_y, self._bound_x

    def drawPolitical(self, scale_len=0):
        self._map.drawcountries(linewidth=1.0)
        self._map.drawcoastlines(linewidth=1.0)

        if not hasattr(self._map, 'counties'):
            self._map.readshapefile("countyp020", 'counties', linewidth=0.5)
        else:
            for county, data in izip(self._map.counties_info, self._map.counties):
                if not self._state_list or county['STATE'] in self._state_list:
                    pylab.gca().add_patch(Polygon(data, ec='k', fc='none', linewidth=0.5))

        if not hasattr(self._map, 'states'):
            self._map.readshapefile("tl_2011_us_state", 'states', linewidth=1.5)
        else:
            for state, data in izip(self._map.states_info, self._map.states):
                if not self._state_list or state['STUSPS'] in self._state_list:
                    pylab.gca().add_patch(Polygon(data, ec='k', fc='none', linewidth=1.5))

        if scale_len > 0:
            pylab.plot([scale_len * 1000, scale_len * 2000], [scale_len * 1000, scale_len * 1000], 'k-', lw=2.5)
            pylab.text(scale_len * 1500, scale_len * 800, "%d km" % scale_len, ha='center', va='top')
        return

    def __call__(self, *args, **kwargs):
        if 'inverse' in kwargs and kwargs['inverse']:
            return self._map(*[ a + d for a, d in zip(args, [1.5 * self._dx, 1.5 * self._dy]) ], **kwargs)
        else:
            return [ c - d for c, d in zip(self._map(*args, **kwargs), [1.5 * self._dx, 1.5 * self._dy]) ]

goshen_1km_grid = lambda bounds=None: \
    ExperimentGrid(projection='lcc', resolution='l',
    grid_size=(255, 255), grid_spacing=(1000, 1000),
    ctr_lat=41.61975, ctr_lon=-104.34843, 
    lat_0=40.619985, lon_0=-107.344, lat_1=30., lat_2=60.,
    state_list=['CO', 'NE', 'WY'], bounds=bounds)

goshen_3km_grid = lambda bounds=None: \
    ExperimentGrid(projection='lcc', resolution='l',
    grid_size=(411, 411), grid_spacing=(3000, 3000),
    lat_0=40.619985, lon_0=-107.344, lat_1=30., lat_2=60.,
    state_list=['CO', 'NE', 'WY', 'SD', 'MT', 'ID', 'UT', 'NV', 'AZ', 'NM', 'TX', 'OK', 'KS' ], bounds=bounds)

if __name__ == "__main__":
    import matplotlib
    matplotlib.use('agg')
    _imports()

    grid = goshen_1km_grid()
    grid.drawPolitical()
    pylab.savefig("grid_test.png")
