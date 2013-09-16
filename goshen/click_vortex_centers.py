
#from util import loadAndInterpolateEnsemble, goshen_1km_proj, goshen_1km_gs, setupMapProjection
from dataload import loadEnsemble
from grid import goshen_1km_grid
from temporal import goshen_1km_temporal
from computeQuantities import computeVorticity, theta2Temperature

import matplotlib.pyplot as pylab
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon

import numpy as np

from itertools import izip
import glob
from datetime import datetime, timedelta
import os
import cPickle

def sliceShape(slice_obj, shape):
    new_slice = []

    def getBounds(b, length):
        if b.start is None: start_x = 0
        else: start_x = b.start

        if b.stop is None: end_x = length
        else: end_x = b.stop

        if b.step is None: step_x = 1
        else: step_x = b.step

        return (end_x - start_x) / step_x

    for sl, sh in zip(slice_obj, shape):
        new_slice.append(getBounds(sl, sh))
    return tuple(new_slice)

class ClickMaxState(object):
    def __init__(self, exp_name, vort_data, ens_members, ens_times, grid):
        self._exp_name = exp_name
        self._current_ens_member = 0
        self._current_time_step = 0
        self._click_dict = {}

        self._loadClicks()

        self._grid = grid
        self._ens_members = ens_members
        if type(self._ens_members) in [ int ]:
            self._ens_members = [ i + 1 for i in range(self._ens_members) ]
        self._ens_times = ens_times

#       gs_x, gs_y = goshen_1km_gs
#       lbx, lby = bounds
#       self._lbx = lbx.start * gs_x
#       self._lby = lby.start * gs_y

        self._setCountyImage()

#       vort_shape = list(vort_data.shape[:2])
#       if bounds is None:
#           vort_shape.extend(vort_data.shape[2:])
#       else:
#           bound_x, bound_y = sliceShape(bounds, vort_data.shape[2:])
#           vort_shape.extend([bound_y, bound_x])

        vort_index = [ slice(None), slice(None) ]
        vort_index.extend(self._grid.getBounds()[::-1])
        
        vort_shape = sliceShape(vort_index, vort_data.shape)

        self._vort_data = np.empty(vort_shape, dtype=vort_data.dtype)
        for var in vort_data.dtype.fields.iterkeys():
            self._vort_data[var] = vort_data[var][vort_index]

        fig = pylab.figure(1)
        self._colorbar = None
        self._vortex_point = None
        cid = fig.canvas.mpl_connect('button_press_event', self.click)
        self.plotMap()

        pylab.show()

        self._saveClicks()
        return

    def click(self, event):
        if event.button == 1:
            self._advanceImage()
        elif event.button == 2:
            ens_tuple = self.__buildExpTuple()
            if event.xdata is not None and event.ydata is not None:
                lby, lbx = [ bnd.start for bnd in self._grid.getBounds() ]
                self._click_dict[ens_tuple] = (event.xdata + lbx, event.ydata + lby)
                print "Placing vortex for", ens_tuple, "at", self._click_dict[ens_tuple], "..."
            else:
                if ens_tuple in self._click_dict:
                    del self._click_dict[ens_tuple]
                    print "Deleting vortex for", ens_tuple
            self._advanceImage()
        elif event.button == 3:
            self._regressImage()
        return

    def plotMap(self, stride=1):
#       pylab.clf()
        pylab.cla()

        vort = self._vort_data['vort'][self._current_ens_member, self._current_time_step]
        u = self._vort_data['u'][self._current_ens_member, self._current_time_step]
        v = self._vort_data['v'][self._current_ens_member, self._current_time_step]

        mag_wind = np.hypot(u, v)
        u_norm = u / mag_wind
        v_norm = v / mag_wind

        thin = (slice(None, None, stride), slice(None, None, stride))

        xs, ys = self._grid.getXY() 

        cs = pylab.contourf(xs, ys, vort, levels=np.arange(-0.015, 0.0175, 0.0025))
        if self._colorbar is not None:
            pylab.gcf().delaxes(pylab.gcf().axes[1])
            pylab.gcf().subplots_adjust(right=0.90)
        self._colorbar = pylab.colorbar()

        pylab.quiver(xs[thin], ys[thin], u[thin], v[thin])

#       self._grid.drawstates(linewidth=1.5)
#       self._grid.drawcountries(linewidth=1.5)
#       self._grid.drawcoastlines(linewidth=1.5)

        pylab.gcf().canvas.restore_region(self._county_image)
        pylab.gcf().canvas.blit(pylab.gca().bbox)

#       if self._vortex_point is not None:
#           pylab.gca().lines.remove(self._vortex_point)
#           self._vortex_point = None

        ens_tuple = self.__buildExpTuple()
        lby, lbx = [ bnd.start for bnd in self._grid.getBounds() ]

        if ens_tuple in self._click_dict:
            x, y = self._click_dict[ens_tuple]
            self._vortex_point = pylab.plot(x - lbx, y - lby, 'ko')[0]

        pylab.title("Member %s, Time %s" % (self._ens_members[self._current_ens_member], self._ens_times[self._current_time_step]), backgroundcolor='w')

        begin = datetime.utcnow()
        pylab.gcf().canvas.draw()
        print "Time to flush canvas:", datetime.utcnow() - begin
        return

    def _setCountyImage(self):

#       if not hasattr(map, 'counties'):
#           self._map.readshapefile("countyp020", 'counties', linewidth=0.5)
#       else:
#           for county, data in izip(map.counties_info, map.counties):
#               if county['STATE'] in ['NE', 'WY', 'CO']:
#                   pylab.gca().add_patch(Polygon(data, ec='k', fc='none', linewidth=0.5))

        self._grid.drawPolitical()

        pylab.gcf().patch.set_alpha(0.0)
        pylab.gca().patch.set_alpha(0.0)

        self._county_image = pylab.gcf().canvas.copy_from_bbox(pylab.gca().bbox)

        pylab.gcf().patch.set_alpha(1.0)
        pylab.gca().patch.set_alpha(1.0)
        return

    def _advanceImage(self):
        if self._current_ens_member == self._vort_data.shape[0] - 1 and self._current_time_step == self._vort_data.shape[1] - 1:
            pass
        else:
            if self._current_time_step + 1 == self._vort_data.shape[1]:
                self._current_time_step = 0
                self._current_ens_member += 1
            else:
                self._current_time_step += 1
            self.plotMap()
        return

    def _regressImage(self):
        if self._current_ens_member == 0 and self._current_time_step == 0:
            pass
        else:
            if self._current_time_step - 1 < 0:
                self._current_time_step = self._vort_data.shape[1] - 1
                self._current_ens_member -= 1
            else:
                self._current_time_step -= 1
            self.plotMap()
        return

    def _saveClicks(self):
        cPickle.dump(self._click_dict, open(self.__buildClickFilename(), 'w'), -1)
        return

    def _loadClicks(self):
        file_name = self.__buildClickFilename()
        if os.path.exists(file_name):
            self._click_dict = cPickle.load(open(file_name, 'r'))
        return

    def __buildClickFilename(self):
        return "vortex_centers_%s.pkl" % self._exp_name

    def __buildExpTuple(self):
        ens = self._ens_members[self._current_ens_member]
        time = self._ens_times[self._current_time_step]
        return (ens, time)

def getVorticity(**kwargs):
    vort_data = np.empty(kwargs['u'].shape, dtype=[('u', np.float32), ('v', np.float32), ('vort', np.float32)])

    vort_data['u'] = kwargs['u']
    vort_data['v'] = kwargs['v']
    vort_data['vort'] = computeVorticity(**kwargs)
    return vort_data

def main():
    bounds = (slice(100, 180), slice(90, 170))
    exp_base = "/caps2/tsupinie/"
    exp_name = "1kmf-r0h=6km-bc7dBZ,5ms"
    height = 2000

    exp_tag = "-".join(exp_name.split("-")[1:])
    n_ens_members = 40
    grid = goshen_1km_grid(bounds=bounds)
    temp = goshen_1km_temporal(start=14400)

#   proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs, bounds)
#   map = Basemap(**proj)

#   fcst_files = glob.glob("%s/%s/ena???.hdf014400" % (exp_base, exp_name))
#   fcst_files.extend(glob.glob("%s/1km-control-%s/ena???.hdf01[5678]*" % (exp_base, exp_name)))

#   vort, ens_members, times = loadAndInterpolateEnsemble(fcst_files, ['u', 'v', 'dx', 'dy'], getVorticity, "%s/%s/ena001.hdfgrdbas" % exp_base, exp_name, { 'z':height })
    vort = loadEnsemble("%s/%s" % (exp_base, exp_name), n_ens_members, temp.getTimes(), (['u', 'v', 'p', 'pt', 'z', 'dx', 'dy'], getVorticity), { 'z':height }, agl=True)

    cms = ClickMaxState("%s-%dm" % (exp_tag, height), vort, n_ens_members, temp, grid)

#   pylab.plot(np.random.random_sample(10), np.random.random_sample(10), 'ko')

#   pylab.xlim(0, 1)
#   pylab.ylim(0, 1)

    return

if __name__ == "__main__":
    main()
