
import matplotlib
matplotlib.use('agg')
import pylab
import numpy as np

import Nio as nio

import cPickle

flux_boxes = {
    'mod-05XP':[
        (slice(119, 140), slice(111, 132)),
        (slice(119, 140), slice(116, 137)),
        (slice(117, 138), slice(120, 141)),
        (slice(117, 138), slice(124, 145)),
        (slice(116, 137), slice(127, 148)),
        (slice(116, 137), slice(131, 152)),
        (slice(115, 136), slice(134, 155)),
        (slice(114, 135), slice(137, 158)),
        (slice(113, 134), slice(140, 161)),
        (slice(112, 133), slice(143, 164)),
        (slice(111, 132), slice(146, 167)),
        (slice(110, 131), slice(149, 170)),
        (slice(109, 130), slice(152, 173)),
    ],

    'mm':[
        (slice(116, 137), slice(112, 133)),
        (slice(117, 138), slice(115, 136)),
        (slice(117, 138), slice(119, 140)),
        (slice(117, 138), slice(122, 143)),
        (slice(116, 137), slice(125, 146)),
        (slice(115, 136), slice(129, 150)),
        (slice(113, 134), slice(131, 152)),
        (slice(112, 133), slice(133, 154)),
        (slice(111, 132), slice(135, 156)),
        (slice(110, 131), slice(137, 158)),
        (slice(109, 130), slice(139, 160)),
        (slice(108, 129), slice(141, 162)),
        (slice(107, 128), slice(143, 164)),
    ],

    'no-mm':[
        (slice(118, 139), slice(110, 131)),
        (slice(117, 138), slice(115, 136)),
        (slice(116, 137), slice(118, 139)),
        (slice(116, 137), slice(121, 142)),
        (slice(115, 136), slice(124, 145)),
        (slice(114, 135), slice(128, 149)),
        (slice(111, 132), slice(129, 150)),
        (slice(110, 131), slice(131, 152)),
        (slice(109, 130), slice(133, 154)),
        (slice(108, 129), slice(135, 156)),
        (slice(107, 128), slice(137, 158)),
        (slice(106, 127), slice(139, 160)),
        (slice(105, 126), slice(141, 162)),
    ]
}

def plotMassFluxProfile(profiles, z_coords, names, title, file_name):
    pylab.figure()

    for prof, zs, name in zip(profiles, z_coords, names):
        pylab.plot(prof, zs / 1000., label=name)

    pylab.xlabel(r"Updraft Mass Flux ($10^6$ kg s$^{-1}$)", size='large')
    pylab.xlim(0, 1.2e9)
    y_lb, y_ub = pylab.ylim()
    pylab.ylim([y_lb, 15])
    locs, labels = pylab.xticks()
    pylab.xticks(locs, (locs / 1e6).astype(int), size='large')

    pylab.ylabel("Height (km)", size='large')
    pylab.yticks(size='large')

    pylab.legend()
    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def main():
    exp_names = {'mod-05XP':"MM + MWR05XP", 'mm':"MM", 'no-mm':"No MM"}
    exp_order = [ 'mod-05XP', 'mm', 'no-mm' ]

    times = np.arange(14400, 18300, 300)
    mass_flux = []
    for exp in exp_order:
        mass_flux.append(cPickle.load(open("mass_flux_%s.pkl" % exp, 'r')))

    grdbas = nio.open_file("/caps1/tsupinie/1km-control-mod-05XP/ena001.hdfgrdbas", mode='r', format='hdf')
    zp = grdbas.variables['zp'][:]

    for wdt, time in enumerate(times):
        column_y = [ (flux_boxes[e][wdt][0].start + flux_boxes[e][wdt][0].stop) / 2 for e in exp_order ]
        column_x = [ (flux_boxes[e][wdt][1].start + flux_boxes[e][wdt][1].stop) / 2 for e in exp_order ]

        z_coords = [ zp[:, y, x] for y, x in zip(column_y, column_x) ]

        plotMassFluxProfile([ m[wdt] for m in mass_flux ], z_coords, [ exp_names[e] for e in exp_order ], "Mass Flux profiles at $t$ = %d" % time, "mass_flux_%06d.png" % time)
    return

if __name__ == "__main__":
    main()
