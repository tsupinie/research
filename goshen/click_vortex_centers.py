
from util import loadAndInterpolateEnsemble, goshen_1km_proj, goshen_1km_gs, setupMapProjection
from computeQuantities import computeVorticity

import matplotlib.pyplot as pylab
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon

import numpy as np

from itertools import izip
import glob
from datetime import datetime, timedelta
import os
import cPickle

class ClickMaxState(object):
    def __init__(self, exp_name, vort_data, ens_members, ens_times, map, bounds=None):
        self._exp_name = exp_name
        self._current_ens_member = 0
        self._current_time_step = 0
        self._click_dict = {}

        self._loadClicks()

        self._map = map
        self._ens_members = ens_members
        self._ens_times = ens_times

        gs_x, gs_y = goshen_1km_gs
        lbx, lby = bounds
        self._lbx = lbx.start * gs_x
        self._lby = lby.start * gs_y

        self._setCountyImage()

        vort_shape = list(vort_data.shape[:2])
        if bounds is None:
            vort_shape.extend(vort_data.shape[2:])
        else:
            def getBounds(b, length):
                if b.start is None: start_x = 0
                else: start_x = b.start

                if b.stop is None: end_x = length
                else: end_x = b.stop

                return end_x - start_x

            bound_y = getBounds(bounds[1], vort_data.shape[2])
            bound_x = getBounds(bounds[0], vort_data.shape[3])

            vort_shape.extend([bound_y, bound_x])

        vort_index = [ slice(None), slice(None) ]
        vort_index.extend(bounds[::-1])
        
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
                self._click_dict[ens_tuple] = (event.xdata + self._lbx, event.ydata + self._lby)
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

        thin = (slice(None, None, stride), slice(None, None, stride))

        nx, ny = vort.shape
        gs_x, gs_y = goshen_1km_gs
        xs, ys = np.meshgrid(gs_x * np.arange(nx), gs_y * np.arange(ny))

        cs = pylab.contourf(xs, ys, vort)
        if self._colorbar is not None:
            pylab.gcf().delaxes(pylab.gcf().axes[1])
            pylab.gcf().subplots_adjust(right=0.90)
        self._colorbar = pylab.colorbar()

        pylab.quiver(xs[thin], ys[thin], u[thin], v[thin])

        self._map.drawstates(linewidth=1.5)
        self._map.drawcountries(linewidth=1.5)
        self._map.drawcoastlines(linewidth=1.5)

        pylab.gcf().canvas.restore_region(self._county_image)
        pylab.gcf().canvas.blit(pylab.gca().bbox)

#       if self._vortex_point is not None:
#           pylab.gca().lines.remove(self._vortex_point)
#           self._vortex_point = None

        ens_tuple = self.__buildExpTuple()
        if ens_tuple in self._click_dict:
            x, y = self._click_dict[ens_tuple]
            self._vortex_point = pylab.plot(x - self._lbx, y - self._lby, 'ko')[0]

        pylab.title("Member %s, Time %s" % (self._ens_members[self._current_ens_member], self._ens_times[self._current_time_step]), backgroundcolor='w')

        begin = datetime.utcnow()
        pylab.gcf().canvas.draw()
        print "Time to flush canvas:", datetime.utcnow() - begin
        return

    def _setCountyImage(self):

        if not hasattr(map, 'counties'):
            self._map.readshapefile("countyp020", 'counties', linewidth=0.5)
        else:
            for county, data in izip(map.counties_info, map.counties):
                if county['STATE'] in ['NE', 'WY', 'CO']:
                    pylab.gca().add_patch(Polygon(data, ec='k', fc='none', linewidth=0.5))

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
            if self._current_ens_member + 1 == self._vort_data.shape[0]:
                self._current_ens_member = 0
                self._current_time_step += 1
            else:
                self._current_ens_member += 1
            self.plotMap()
        return

    def _regressImage(self):
        if self._current_ens_member == 0 and self._current_time_step == 0:
            pass
        else:
            if self._current_ens_member - 1 < 0:
                self._current_ens_member = self._vort_data.shape[0] - 1
                self._current_time_step -= 1
            else:
                self._current_ens_member -= 1
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
    exp_base = "/caps1/tsupinie/"
    exp_name = "mod-05XP"
    height = 500

    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs, bounds)
    map = Basemap(**proj)

    fcst_files = glob.glob("%s/1km-control-%s/ena???.hdf014400" % (exp_base, exp_name))
#   fcst_files.extend(glob.glob("%s/1km-control-%s/ena???.hdf01[5678]*" % (exp_base, exp_name)))

    vort, ens_members, times = loadAndInterpolateEnsemble(fcst_files, ['u', 'v', 'dx', 'dy'], getVorticity, "%s/1km-control-20120712/ena001.hdfgrdbas" % exp_base, { 'z':height })

    cms = ClickMaxState("%s-%dm" % (exp_name, height), vort, ens_members, times, map, bounds)

#   pylab.plot(np.random.random_sample(10), np.random.random_sample(10), 'ko')

#   pylab.xlim(0, 1)
#   pylab.ylim(0, 1)

    return

if __name__ == "__main__":
    main()
