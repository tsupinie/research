
from util import loadAndInterpolateEnsemble, setupMapProjection, goshen_1km_proj, goshen_1km_gs, drawPolitical
from computeQuantities import computeVorticity

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap

import numpy as np

import glob
from datetime import datetime, timedelta

def vorticityWinds(**kwargs):
    vortWinds = np.empty(kwargs['u'].shape, dtype=[('u', np.float32), ('v', np.float32), ('vort', np.float32)])

    vortWinds['u'] = kwargs['u']
    vortWinds['v'] = kwargs['v']
    vortWinds['vort'] = computeVorticity(**kwargs)
    return vortWinds

def plotVorticity(vorticity, winds, map, title, file_name, stride=4):
    u, v = winds
    pylab.figure()

    thin_data = tuple([ slice(None, None, stride) ] * 2)
    gs_x, gs_y = goshen_1km_gs
    nx, ny = vorticity.shape
    xs, ys = np.meshgrid(gs_x * np.arange(nx), gs_y * np.arange(ny))

    pylab.contourf(xs, ys, vorticity, levels=np.arange(0.0075, 0.034, 0.0025))
    pylab.colorbar()
    pylab.quiver(xs[thin_data], ys[thin_data], u[thin_data], v[thin_data])

    drawPolitical(map)

    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def main():
    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs)
    map = Basemap(**proj)
    exp_order = [ 'no-mm', 'mm', 'mod-05XP' ]
    base_path = "/caps1/tsupinie/1km-control-emeanf/"

    files = glob.glob("%s/ena???.hdf0*" % base_path)
    
    vort_all, ens_members, times = loadAndInterpolateEnsemble(files, [ 'u', 'v', 'dx', 'dy' ], vorticityWinds, "%s/ena001.hdfgrdbas" % base_path, { 'z':2000 })

    print vort_all.shape

    for lde, ens in enumerate(ens_members):
        for wdt, time in enumerate(times):
            plotVorticity(vort_all['vort'][lde, wdt], (vort_all['u'][lde, wdt], vort_all['v'][lde, wdt]), map, "Vorticity", "det_vort_%s_%06d.png" % (exp_order[lde], time))
    return

if __name__ == "__main__":
    main()
