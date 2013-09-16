
from legacy import goshen_1km_proj, goshen_1km_gs, setupMapProjection

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap

import cPickle

def reorganizeCenters(centers, times):
    lines_x = []
    lines_y = []

    def uniquify(L):
        seen = {}
        unique = []

        for item in L:
            if item not in seen:
                unique.append(item)
                seen[item] = True

        return unique

    all_members, all_times = zip(*centers.keys())

    ens_members = sorted(uniquify(all_members))

    for ens_str in uniquify(ens_members):
        line_x = np.empty((len(times), ))
        line_y = np.empty((len(times), ))

        for wdt, t_ens in enumerate(times):
            if (ens_str, t_ens) in centers:
                x, y = centers[(ens_str, t_ens)]
                line_x[wdt] = x
                line_y[wdt] = y
            else:
                line_x[wdt] = np.nan
                line_y[wdt] = np.nan

        lines_x.append(line_x)
        lines_y.append(line_y)

    return lines_x, lines_y

def reorganizeCentersDumb(centers, n_ens_members, times):
    lines_x = []
    lines_y = []

    for n_ens in xrange(n_ens_members):
        ens_str = "%03d" % (n_ens + 1)

        line_x = np.empty(times.shape)
        line_y = np.empty(times.shape)

        for wdt, t_ens in enumerate(times):
            if (ens_str, t_ens) in centers:
                x, y = centers[(ens_str, t_ens)]
                line_x[wdt] = x
                line_y[wdt] = y
            else:
                line_x[wdt] = np.nan
                line_y[wdt] = np.nan

        lines_x.append(line_x)
        lines_y.append(line_y)
    return lines_x, lines_y

def main():
    exp_name = "mm"
    height = 1000
    n_ens_members = 40
    times = np.arange(14400, 18300, 300)

    centers = cPickle.load(open("vortex_centers_%s-%dm.pkl" % (exp_name, height), 'r'))
    lines_x, lines_y = reorganizeCenters(centers, times) #, n_ens_members, times)
#   lines_x, lines_y = reorganizeCentersDumb(centers, n_ens_members, times)

    bounds = (slice(100, 180), slice(90, 170))
    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs, bounds)
    map = Basemap(**proj)

    gs_x, gs_y = goshen_1km_gs
    lb_x = bounds[0].start * gs_x
    lb_y = bounds[1].start * gs_y

    tornado_track_x, tornado_track_y = map(*zip(*((41.63,-104.383), (41.6134,-104.224)))[::-1])
    pylab.plot(tornado_track_x, tornado_track_y, 'r-')

    for line_x, line_y in zip(lines_x, lines_y):
        pylab.plot(line_x - lb_x, line_y - lb_y, 'ko-', markersize=2, linewidth=1)

    map.drawstates(linewidth=1.5)
    map.drawcountries(linewidth=1.5)
    map.drawcoastlines(linewidth=1.5)
    map.readshapefile("countyp020", 'counties', linewidth=0.5)

    pylab.savefig("vortex_center_swath_%s-%dm.png" % (exp_name, height))
    return

if __name__ == "__main__":
    main()
