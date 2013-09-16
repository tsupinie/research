
import cPickle, glob
from operator import mul

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib.ticker import FixedFormatter, FixedLocator
from matplotlib.colors import LinearSegmentedColormap

from grid import goshen_1km_grid
from temporal import goshen_1km_temporal

confusion_cmap_dict = {
    'red':(
        (0.0, 1.0, 1.0), 
        (0.2, 1.0, 0.5), # Missing to Correct Negative
        (0.4, 0.5, 1.0), # Correct Negative to False Alarm
        (0.6, 1.0, 1.0), # False Alarm to Miss
        (0.8, 1.0, 0.0), # Miss to Hit
        (1.0, 0.0, 0.0),
    ),
    'green':(
        (0.0, 1.0, 1.0),
        (0.2, 1.0, 0.5), 
        (0.4, 0.5, 0.0),
        (0.6, 0.0, 1.0),
        (0.8, 1.0, 0.5),
        (1.0, 0.5, 0.5),
    ),
    'blue':(
        (0.0, 1.0, 1.0),
        (0.2, 1.0, 0.5), 
        (0.4, 0.5, 1.0),
        (0.6, 1.0, 0.0),
        (0.8, 0.0, 0.0),
        (1.0, 0.0, 0.0),
    ),
}
confusion_cmap = LinearSegmentedColormap('confusion', confusion_cmap_dict, 256)

def plotConfusion(confusion, grid, title, file_name, inset=None, fudge=16):
    pylab.figure()
    axmain = pylab.axes()

    tick_labels = [ "Missing", "Correct\nNegative", "False\nAlarm", "Miss", "Hit" ]
    min_label = -1

    xs, ys = grid.getXY()

    pylab.pcolormesh(xs, ys, confusion, cmap=confusion_cmap, vmin=min_label, vmax=(min_label + len(tick_labels) - 1))

    tick_locs = np.linspace(-1, min_label + len(tick_labels) - 2, len(tick_labels))
    tick_locs += (tick_locs[1] - tick_locs[0]) / 2
    bar = pylab.colorbar()
    bar.locator = FixedLocator(tick_locs)
    bar.formatter = FixedFormatter(tick_labels)
    pylab.setp(pylab.getp(bar.ax, 'ymajorticklabels'), fontsize='large')
    bar.update_ticks()

    grid.drawPolitical()

    if inset:
        lb_y, lb_x = [ b.start for b in inset ]
        ub_y, ub_x = [ b.stop + fudge for b in inset ]

        inset_exp = (slice(lb_y, ub_y), slice(lb_x, ub_x))

        axins = zoomed_inset_axes(pylab.gca(), 2, loc=4)
        pylab.sca(axins)

        pylab.pcolormesh(xs[inset_exp], ys[inset_exp], confusion[inset_exp], cmap=confusion_cmap, vmin=min_label, vmax=(min_label + len(tick_labels) - 1))
        grid.drawPolitical()

        gs_x, gs_y = grid.getGridSpacing()

        pylab.xlim([lb_x * gs_x, (ub_x - 1) * gs_x])
        pylab.ylim([lb_y * gs_y, (ub_y - 1) * gs_y])

        mark_inset(axmain, axins, loc1=1, loc2=3, fc='none', ec='k')

    pylab.sca(axmain)
    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def plotGraph(xs, ys, title, file_name):
    pylab.figure()
    for key, y in ys.iteritems():
        pylab.plot(xs.getTimes(), y, label=key)

    pylab.legend(loc=2)
    pylab.xticks(xs.getTimes(), xs.getStrings("%H%M", aslist=True), rotation=30)
    pylab.xlim(xs.getTimes()[0], xs.getTimes()[-1])
    pylab.ylim(0., 0.24)
    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def findSize(bbox):
    return reduce(mul, [ b.stop - b.start for b in bbox ])

def main():
    refl_threshold = 20
    vel_threshold = 20

    confusion_refl = cPickle.load(open("all_confusion_1km_%ddBZ.pkl" % (refl_threshold)))
    confusion_vel = cPickle.load(open("all_confusion_1km_%dms.pkl" % (vel_threshold)))

    missing_data = -1
    fa_value = 1
    miss_value = 2

    grid = goshen_1km_grid()
    temp = goshen_1km_temporal(start=14400)

    bbox_files = glob.glob("bbox*.pkl")
    bbox_buffer = 10
    bbox_offsets = [0, 10, 0]
    bboxes = {}
    for bfile in bbox_files:
        root, ext = bfile.split(".")
        bits = root.split("_")
        key = "1kmf-%s" % bits[-1]
        bboxes[key] = cPickle.load(open(bfile, 'r'))

        bboxes[key] = (slice(0, 14),) + tuple( slice(b.start - bbox_buffer + o, b.stop + bbox_buffer + o) for b, o in zip(bboxes[key][1:], bbox_offsets[1:]) )

    misses       = dict( (e, [ (c[bboxes[e]] == miss_value).sum() / float((c[bboxes[e]] != missing_data).sum()) for c in r['KCYS'] ]) for e, r in confusion_refl.iteritems() )
    false_alarms = dict( (e, [ (c[bboxes[e]] == fa_value).sum() / float((c[bboxes[e]] != missing_data).sum()) for c in r['KCYS'] ]) for e, r in confusion_refl.iteritems() )

    plotGraph(temp, misses, "Misses", "all_misses_%ddBZ_thd0.5.png" % refl_threshold)
    plotGraph(temp, false_alarms, "False Alarms", "all_false_alarms_%ddBZ_thd0.5.png" % refl_threshold)

    for exp, r_confusion in confusion_refl.iteritems():
        for confusion, time_sec, time_str in zip(r_confusion['KCYS'], temp.getTimes(), temp.getStrings('%H%M')):
            plotConfusion(confusion[0], grid, "Contingency Plot for KCYS at %s" % time_str, "confusion_KCYS_%s_tilt00_%06d.png" % (exp, time_sec), inset=bboxes[exp][1:])

    return

if __name__ == "__main__":
    main()
