
from util import setupMapProjection, goshen_3km_proj, goshen_3km_gs, loadAndInterpolateEnsemble, drawPolitical
from computeQuantities import computePMSL, theta2Temperature, qv2Dewpoint

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap

import glob

def plotSurface(pt, td, winds, map, stride, title, file_name):
    pylab.figure()
    pylab.axes((0.05, 0.025, 0.9, 0.9))

    u, v = winds

    nx, ny = pt.shape
    gs_x, gs_y = goshen_3km_gs
    xs, ys = np.meshgrid(gs_x * np.arange(nx), gs_y * np.arange(ny))   

    data_thin = tuple([ slice(None, None, stride) ] * 2)

    td_cmap = matplotlib.cm.get_cmap('Greens')
    td_cmap.set_under('#ffffff')
    pylab.contourf(xs, ys, td, levels=np.arange(40, 80, 5), cmap=td_cmap)
    pylab.colorbar()
    CS = pylab.contour(xs, ys, pt, colors='r', linestyles='-', linewidths=1.5, levels=np.arange(288, 324, 4))
    pylab.clabel(CS, inline_spacing=0, fmt="%d K", fontsize='x-small')
    pylab.quiver(xs[data_thin], ys[data_thin], u[data_thin], v[data_thin])

    drawPolitical(map, scale_len=75)

    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def getSurface(**kwargs):
    sfc = np.empty(kwargs['u'].shape, dtype=[('u', np.float32), ('v', np.float32), ('pt', np.float32), ('td', np.float32)])

    sfc['u'] = kwargs['u']
    sfc['v'] = kwargs['v']
    sfc['pt'] = kwargs['pt'] 
    sfc['td'] = 9. / 5. * (qv2Dewpoint(**kwargs) - 273.15) + 32
    return sfc

def main():
    base_path = "/caps1/tsupinie/3km-ensemble-20120712/"
    files = glob.glob("%s/ena???.hdf010800" % base_path)

    proj = setupMapProjection(goshen_3km_proj, goshen_3km_gs)
    map = Basemap(**proj)

    sfc, ens_members, ens_times = loadAndInterpolateEnsemble(files, ['p', 'pt', 'qv', 'u', 'v'], getSurface, "3km-orig-ensemble/ena011.hdfgrdbas",
        {'z':10}, wrap=True, agl=True)

    u = sfc['u'].mean(axis=0)[0]
    v = sfc['v'].mean(axis=0)[0]
    pt = sfc['pt'].mean(axis=0)[0]
    td = sfc['td'].mean(axis=0)[0]

    plotSurface(pt, td, (u, v), map, 8, "Surface plot at 2100 UTC", "sfc_3km.png")

if __name__ == "__main__":
    main()
