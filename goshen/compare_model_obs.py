
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

from grid import goshen_3km_grid
from arpsmodelobs import ARPSModelObsFile
from util import coneHeight

import cPickle

def main():
    arps_mo = ARPSModelObsFile("/caps1/tsupinie/3km-fixed-radar/KCYSan014400", (2, 12))
    python_mo = cPickle.load(open("3km-fixed-radar/eomean.pkl014400", 'r'))

    grid = goshen_3km_grid()#bounds=(slice(242, 327), slice(199, 284)))
    radar_x, radar_y = grid(arps_mo._radar_lon, arps_mo._radar_lat)
    print radar_x, arps_mo._radar_x
    print radar_y, arps_mo._radar_y
    xs, ys = grid.getXY()
    bounds = grid.getBounds()

    Re = 6371000.
    range = np.hypot(xs - radar_x, ys - radar_y)
    slant_range = Re * np.tan(range / Re)
    range_mask = (range < arps_mo._range_max) & (range > arps_mo._range_min)

    height_error = arps_mo._height[0][bounds] - coneHeight(range, arps_mo._elev_angles[0] + 0.06, 1867.)
    height_error_masked = np.where(range < arps_mo._range_max - 12000, height_error, np.nan)
    print np.nanmin(height_error_masked), np.nanmax(height_error_masked)

    pylab.figure(figsize=(12,6))

    pylab.subplot(121)
    pylab.pcolormesh(xs, ys, arps_mo['Z'][0][bounds], vmin=10, vmax=80, cmap=pylab.jet())
#   pylab.pcolormesh(xs, ys, height_error_masked, vmin=-200, vmax=200, cmap=matplotlib.cm.get_cmap('RdBu'))
    pylab.colorbar()
    pylab.title("ARPS")
    grid.drawPolitical()

    pylab.subplot(122)
    pylab.pcolormesh(xs, ys, python_mo['Z'][0, 0][bounds], vmin=10, vmax=80, cmap=matplotlib.cm.get_cmap('jet'))
    pylab.title("Python")
    grid.drawPolitical()

    pylab.savefig("model_obs_comparison.png")
    pylab.close()

    cPickle.dump(arps_mo._height, open("beam_height.pkl", 'w'), -1)

    return

if __name__ == "__main__":
    main()
