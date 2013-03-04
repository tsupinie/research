
from util import loadAndInterpolateEnsemble, setupMapProjection, goshen_1km_proj, goshen_1km_gs, drawPolitical, flux_boxes
from computeQuantities import computeReflectivity

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

import glob, argparse

def makeplot(refl, winds, w, stride, map, gs, title, file_name, box=None):
    pylab.figure()
    axmain = pylab.axes((0, 0.025, 1, 0.9))

    gs_x, gs_y = gs
    nx, ny = refl.shape
    xs, ys = np.meshgrid(gs_x * np.arange(nx), gs_y * np.arange(ny))

    pylab.contourf(xs, ys, refl, levels=np.arange(10, 80, 10))
    pylab.colorbar()
    pylab.contour(xs, ys, w, levels=np.arange(-10, 0, 2), colors='#666666', style='--')
    pylab.contour(xs, ys, w, levels=np.arange(2, 12, 2), colors='#666666', style='-')
    
    u, v = winds
    wind_slice = tuple([ slice(None, None, stride) ] * 2)
    pylab.quiver(xs[wind_slice], ys[wind_slice], u[wind_slice], v[wind_slice])

    if box:
        lb_y, lb_x = [ b.start for b in box ]
        ub_y, ub_x = [ b.stop for b in box ]

        box_xs = gs_x * np.array([ lb_x, lb_x, ub_x, ub_x, lb_x])
        box_ys = gs_y * np.array([ lb_y, ub_y, ub_y, lb_y, lb_y])

        map.plot(box_xs, box_ys, '#660099')

        axins = zoomed_inset_axes(pylab.gca(), 4, loc=4)
        pylab.sca(axins)

        pylab.contourf(xs[box], ys[box], refl[box], levels=np.arange(10, 80, 10))
        pylab.contour(xs[box], ys[box], w[box], levels=np.arange(-10, 0, 2), colors='#666666', style='--')
        pylab.contour(xs[box], ys[box], w[box], levels=np.arange(2, 12, 2), colors='#666666', style='-')
        pylab.quiver(xs[box], ys[box], u[box], v[box])

        drawPolitical(map)

        pylab.xlim([lb_x * gs_x, ub_x * gs_x - 1])
        pylab.ylim([lb_y * gs_y, ub_y * gs_y - 1])

        mark_inset(axmain, axins, loc1=1, loc2=3, fc='none', ec='k')

    pylab.sca(axmain)
    drawPolitical(map)

    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def getWind(**kwargs):
    wind = np.empty(kwargs['u'].shape, dtype=[('u', np.float32), ('v', np.float32), ('w', np.float32)])

    wind['u'] = kwargs['u']
    wind['v'] = kwargs['v']
    wind['w'] = kwargs['w']
    return wind

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp-name', dest='exp_name', required=True)

    args = ap.parse_args()

    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs)
    map = Basemap(**proj)

    fcst_files = glob.glob("/caps1/tsupinie/1km-control-%s/ena???.hdf0*" % args.exp_name)

    refl_ens_mean, ens_refl, ens_members, ens_times = loadAndInterpolateEnsemble(fcst_files, ['pt', 'p', 'qr', 'qs', 'qh'], computeReflectivity, "/caps1/tsupinie/1km-control-%s/ena001.hdfgrdbas" % args.exp_name, 
        {'z':1000}, agl=True, aggregator=lambda x: np.mean(x, axis=0))

    ens_winds, ens_members, ens_times = loadAndInterpolateEnsemble(fcst_files, ['u', 'v', 'w'], getWind, "/caps1/tsupinie/1km-control-%s/ena001.hdfgrdbas" % args.exp_name, 
        {'z':1000}, agl=True)
    
    u_mean = ens_winds['u'].mean(axis=0)
    v_mean = ens_winds['v'].mean(axis=0)
    w_mean = ens_winds['w'].mean(axis=0)

    fcst_start_idx = np.where(ens_times == 14400)[0][0]

    for wdt, t_ens in enumerate(ens_times):
        if wdt < fcst_start_idx:
            makeplot(refl_ens_mean[wdt], (u_mean[wdt], v_mean[wdt]), w_mean[wdt], 5, map, goshen_1km_gs, r"Ensemble Mean Reflectivity, $w$, and Wind at $z$ = 1000 m and $t$ = %d s" % t_ens, "images-%s/ens_mean_%06d.png" % (args.exp_name, t_ens))
        else:
            makeplot(refl_ens_mean[wdt], (u_mean[wdt], v_mean[wdt]), w_mean[wdt], 5, map, goshen_1km_gs, r"Ensemble Mean Reflectivity, $w$, and Wind at $z$ = 1000 m and $t$ = %d s" % t_ens, "images-%s/ens_mean_%06d.png" % (args.exp_name, t_ens), 
                box=flux_boxes[args.exp_name][wdt - fcst_start_idx])

    return

if __name__ == "__main__":
    main()
