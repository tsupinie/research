
import numpy as np

from radarobsfile import RadarObsFile
from grid import goshen_1km_grid

import matplotlib
matplotlib.use('agg')
import pylab

def main():
    mwr = RadarObsFile("qc/1km/vols/05XP.20090605.220000")
    wsr = RadarObsFile("qc/1km/KCYS.20090605.220000")

    bounds = (slice(110, 135), slice(118, 143))
    grid = goshen_1km_grid(bounds)    
    bounds = grid.getBounds()
    xs, ys = grid.getXY()
    
    pylab.contourf(xs, ys, wsr['vr'][0][bounds], levels=np.arange(-40, 45, 5), cmap=matplotlib.cm.get_cmap('RdBu_r'))
    pylab.colorbar()
    pylab.contour(xs, ys, mwr['vr'][3][bounds], levels=np.arange(-40, 0, 5), linestyle='--', colors='k')
    pylab.contour(xs, ys, mwr['vr'][3][bounds], levels=[0], linestyle='-', linewidths=1.5, colors='k')
    pylab.contour(xs, ys, mwr['vr'][3][bounds], levels=np.arange(5, 45, 5), linestyle='-', colors='k')
    pylab.contour(xs, ys, mwr['Z'][0][bounds], levels=[10], linewidths=2, colors='k')

    grid.drawPolitical()

    pylab.savefig('overlay.png')
    return

if __name__ == "__main__":
    main()
