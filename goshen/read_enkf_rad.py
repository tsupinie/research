
from util import setupMapProjection, goshen_1km_proj, goshen_1km_gs, goshen_3km_proj, goshen_3km_gs, drawPolitical

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap

import struct
from datetime import datetime, timedelta
from operator import mul
from itertools import izip
import glob

from radarobsfile import RadarObsFile

mask_cmap_dict = {
    'red':(
        (0.0,  1.0, 1.0), 
        (1.0,  0.75, 0.75),
    ),
    'green':(
        (0.0,  1.0, 1.0), 
        (1.0,  0.75, 0.75),
    ),
    'blue':(
        (0.0,  1.0, 1.0), 
        (1.0,  0.75, 0.75),
    ),
}
mask_cmap = LinearSegmentedColormap('confusion', mask_cmap_dict, 256)

def plotRadTilt(plot_data, plot_cmaps, plot_titles, grid, file_name, base_ref=None):
    subplot_base = 220
    n_plots = 4

    xs, ys, gs_x, gs_y, map = grid

    pylab.figure(figsize=(12,12))
    pylab.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02, wspace=0.04)

    for plot in range(n_plots):
        pylab.subplot(subplot_base + plot + 1)
        cmap, vmin, vmax = plot_cmaps[plot]
        cmap.set_under("#ffffff", alpha=0.0)

        pylab.pcolormesh(xs - gs_x / 2, ys - gs_y / 2, plot_data[plot] >= -90, vmin=0, vmax=1, cmap=mask_cmap)

        pylab.pcolormesh(xs - gs_x / 2, ys - gs_y / 2, plot_data[plot], vmin=vmin, vmax=vmax, cmap=cmap)
        pylab.colorbar()

        if base_ref is not None:
            pylab.contour(xs, ys, base_ref, levels=[10.], colors='k', lw=0.5)

#       if plot == 0:
#           pylab.contour(xs, ys, refl_88D[:, :, 0], levels=np.arange(10, 80, 10), colors='#808080', lw=0.5)

#       pylab.plot(radar_x, radar_y, 'ko')

        pylab.title(plot_titles[plot])

        drawPolitical(map)

    pylab.savefig(file_name)
    pylab.close()

    return

def main():
#   files = glob.glob("qc/manual/3km/KCYS.20090605.*")
    files = glob.glob("KCYS.20090605.*")
#   erf_88D = RadarObsFile("qc/manual/1km/KCYS.20090605.215744")
#   domain_bounds = (slice(90, 160), slice(80, 150))
    domain_bounds = (slice(None), slice(None))

    radar_location = (41.56150, -104.298996)

    proj = setupMapProjection(goshen_3km_proj, goshen_3km_gs, bounds=tuple(reversed(domain_bounds)))
    map = Basemap(**proj)

    radar_x, radar_y = map(*reversed(radar_location))

    gs_x, gs_y = goshen_3km_gs

    plot_cmaps = [
        (matplotlib.cm.jet, 10, 80),
        (matplotlib.cm.RdBu, -40, 40),
        (matplotlib.cm.RdYlGn, 0, 10000),
        (matplotlib.cm.RdYlGn, 0, 40000),
    ]

    plot_titles = [
        "Reflectivity (dBZ)",
        "Radial Velocity (m s$^{-1}$)",
        "Height (m)",
        "Range (m)",
    ]

    for file in files:
        print file
        erf = RadarObsFile(file)

        xs, ys = np.meshgrid(gs_x * np.arange(erf._n_gridx), gs_y * np.arange(erf._n_gridy))

        xs = xs[domain_bounds]
        xs -= xs[0, 0]
        ys = ys[domain_bounds]
        ys -= ys[0, 0]

        plot_data = [ 
            erf.reflectivity[domain_bounds],
            erf.radial_velocity[domain_bounds],
            erf.heights[domain_bounds],
            erf.range[domain_bounds],
        ]

        file_name_start = file.rfind("/") + 1
        radar_id = file[file_name_start:(file_name_start + 4)]
        time = file[-6:]

        for ntlt in [ 0 ]: #range(erf._n_tilts):
            plotRadTilt([ p[:, :, ntlt] for p in plot_data ], plot_cmaps, plot_titles, (xs, ys, gs_x, gs_y, map), "%s_enkf_rad_%s.png" % (radar_id, time), base_ref=plot_data[0][:, :, 0])
    return

if __name__ == "__main__":
    main()
