
import cPickle

from util import plot_map, setupMapProjection

import numpy as np

import matplotlib
import pylab

def main():
    grid_spacing = 250
    center_lat, center_lon = 41.61975, -104.34843
    max_runs = cPickle.load(open("250m-det.pkl"))
    keys = sorted(max_runs.keys())

    swath_x_plots, swath_y_plots = 2, 3
    diff_x_plots, diff_y_plots = 3, 4

    ny, nx = max_runs[keys[0]].shape

    proj = {'projection':'lcc', 'resolution':'l',
        'width':(grid_spacing * nx), 'height':(grid_spacing * ny),
        'ctr_lat':center_lat, 'ctr_lon':center_lon, 
        'lat_0':40.61998, 'lon_0':-107.344, 'lat_1':30., 'lon_1':60.
    }

    proj = setupMapProjection(proj, (grid_spacing, grid_spacing))

    cmix_diffs = {}
    for key1 in keys:
        for key2 in keys:
            diff_key = "%s:%s" % (key1, key2)
            diff_key_rev = "%s:%s" % (key2, key1)
            if key1 != key2 and diff_key not in cmix_diffs and diff_key_rev not in cmix_diffs:
                refl1 = np.maximum(max_runs[key1], 10 * np.ones(max_runs[key1].shape))
                refl2 = np.maximum(max_runs[key2], 10 * np.ones(max_runs[key2].shape))
                diff_refl = refl2 - refl1

                cmix_diffs[diff_key] = diff_refl

    pylab.figure(figsize=(16, 10))
    pylab.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02, wspace=0.04)

    subplot_base = swath_x_plots * 100 + swath_y_plots * 10
    for idx, key in enumerate(keys):
        plot_map(max_runs[key], proj, key, color_bar='refl', subplot=(subplot_base + idx + 1), pcolormesh=True)
    pylab.savefig("cmix_swaths.png")

    pylab.figure(figsize=(16, 10))
    pylab.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02, wspace=0.04)

    subplot_base = diff_x_plots * 100 + diff_y_plots * 10
    for idx, key in enumerate(sorted(cmix_diffs.keys())):
        plot_map(cmix_diffs[key], proj, key, color_bar='drefl', subplot=(subplot_base + idx + 1), pcolormesh=True)
    pylab.savefig("cmix_diffs.png")
    return

if __name__ == "__main__":
    main()
