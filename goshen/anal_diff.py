
from legacy import loadAndInterpolateEnsemble, setupMapProjection, decompressVariable, goshen_1km_proj, goshen_1km_gs
from computeQuantities import computeReflectivity

import pylab
from mpl_toolkits.basemap import Basemap

import numpy as np
import Nio as nio

from math import floor, ceil
import glob

forecast_base = "/caps1/tsupinie/1km-control-20120927/"
analysis_base = "/caps1/tsupinie/1km-control-20120927/"
interp_height = 1000

def main(variable, refl_ens_mean):
    radar_elev, radar_lat, radar_lon = 1883, 41.151944, -104.806111
    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs)

    map = Basemap(**proj)
    radar_x, radar_y = map(radar_lon, radar_lat)

    ena_files      = glob.glob("%s/ena???.hdf01[123]*" % analysis_base) #01[123]
    ena_files.extend(glob.glob("%s/ena???.hdf014[14]*" % analysis_base))
    enf_files      = glob.glob("%s/enf???.hdf0*" % forecast_base)

    ena_pt, ens_members, ens_times = loadAndInterpolateEnsemble(ena_files, [variable], lambda **x: x[variable], "/caps1/tsupinie/1km-control-no-ua/ena001.hdfgrdbas", 
        {'z':interp_height}, agl=True, wrap=False)

    enf_pt, ens_members, ens_times = loadAndInterpolateEnsemble(enf_files, [variable], lambda **x: x[variable], "/caps1/tsupinie/1km-control-20120712/ena001.hdfgrdbas", 
        {'z':interp_height}, agl=True, wrap=False)

    ena_mean_pt = ena_pt.mean(axis=0)
    enf_mean_pt = enf_pt.mean(axis=0)

    ena_sprd_pt = ena_pt.std(axis=0, ddof=1)
    enf_sprd_pt = enf_pt.std(axis=0, ddof=1)

    mean_diff = ena_mean_pt - enf_mean_pt
    sprd_diff = ena_sprd_pt - enf_sprd_pt

    dx, dy = goshen_1km_gs
    nx, ny = mean_diff[0].shape

    xs, ys = np.meshgrid(dx * np.arange(nx), dy * np.arange(ny))

    def getLevels(min, max, precision=1):
        levels_max = ceil(float(max) / precision) * precision + precision
        levels_min = floor(float(min) / precision) * precision
        return np.arange(levels_min, levels_max, precision)

    for wdt, ens_t in enumerate(ens_times):
        pylab.figure(figsize=(13, 6))
        print "Minimum/maximum mean difference:", np.nanmin(mean_diff), np.nanmax(mean_diff)
        pylab.subplot(121)
        map.contourf(xs, ys, mean_diff[wdt], levels=getLevels(np.nanmin(mean_diff), np.nanmax(mean_diff), 1))
        pylab.colorbar()
        map.contour(xs, ys, refl_ens_mean[wdt], levels=[20, 40], colors='k')

        map.drawcoastlines(linewidth=1.5)
        map.drawcountries(linewidth=1.5)
        map.drawstates(linewidth=1.0)
        map.readshapefile('countyp020', 'counties', linewidth=0.5)

        print "Minimum/maximum spread difference:", np.nanmin(sprd_diff), np.nanmax(sprd_diff)

        pylab.subplot(122)
        map.contourf(xs, ys, sprd_diff[wdt], levels=getLevels(np.nanmin(sprd_diff), np.nanmax(sprd_diff), 0.05))
        pylab.colorbar()
        map.contour(xs, ys, refl_ens_mean[wdt], levels=[20, 40], colors='k')

        map.drawcoastlines(linewidth=1.5)
        map.drawcountries(linewidth=1.5)
        map.drawstates(linewidth=1.0)
        map.readshapefile('countyp020', 'counties', linewidth=0.5)

        pylab.savefig("mean_sprd_diff_%s_%s_no-rad.png" % (ens_t, variable))
    return

if __name__ == "__main__":
    ena_files      = glob.glob("%s/ena???.hdf01[123]*" % analysis_base)  #01[123]
    ena_files.extend(glob.glob("%s/ena???.hdf014[14]*" % forecast_base))

    refl_ens_mean, ens_refl, ens_members, ens_times = loadAndInterpolateEnsemble(ena_files, ['pt', 'p', 'qr', 'qs', 'qh'], computeReflectivity, "/caps1/tsupinie/1km-control-no-ua/ena001.hdfgrdbas", 
        {'z':interp_height}, agl=True, wrap=False, aggregator=lambda x: np.mean(x, axis=0))

    for var in [ 'u', 'v', 'w', 'pt']:
        main(var, refl_ens_mean)
