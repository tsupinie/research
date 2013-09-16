
import numpy as np

import matplotlib
from matplotlib.ticker import FixedFormatter, FixedLocator
matplotlib.use('agg')
import pylab
from matplotlib.collections import LineCollection
from mpl_toolkits.basemap import Basemap
from shapelib import ShapeFile
#import dbflib

import cPickle
from datetime import datetime, timedelta
import argparse

from legacy import setupMapProjection, goshen_1km_proj, goshen_1km_gs, goshen_3km_proj, goshen_3km_gs, drawPolitical
#from grid import goshen_1km_grid
from plot_vortex_centers import reorganizeCenters

def loadCountyBoundaries(file_name, map):
    shape_file = ShapeFile(file_name)
    n_polys = shape_file.info()[0]

    all_verts_ll = []
    all_lengths = [ 0 ]
    for n_poly in xrange(n_polys):
        shape = shape_file.read_object(n_poly)
        rings = shape.vertices()
        for ring in rings:
            verts = np.array(zip(*ring))
            all_verts_ll.append(verts)
            all_lengths.append(verts.shape[1])

    start = datetime.utcnow()
    vert_glob_ll = np.hstack(tuple(all_verts_ll))
    vert_glob_xy = np.array(zip(*map(*vert_glob_ll)))
    print "Time to map and zip:", datetime.utcnow() - start

    partitions = np.add.accumulate(all_lengths)
    all_verts_xy = []

    for lbound, ubound in zip(partitions[:-1], partitions[1:]):
        all_verts_xy.append( vert_glob_xy[lbound:ubound] )

    return LineCollection(all_verts_xy, antialiaseds=(1,), color='k', lw=0.5)

def removeNones(object, shape):
    obj_t, obj_y, obj_x = object
    if obj_t.start is None: obj_t.start = 0
    if obj_t.stop is None:  obj_t.stop = shape[0]
    if obj_y.start is None: obj_y.start = 0
    if obj_y.stop is None:  obj_y.stop = shape[1]
    if obj_x.start is None: obj_x.start = 0
    if obj_x.stop is None:  obj_x.stop = shape[2]

    return obj_t, obj_y, obj_x

def findSize(object):
    obj_t, obj_y, obj_x = object
    return (obj_t.stop - obj_t.start) * (obj_y.stop - obj_y.start) * (obj_x.stop - obj_x.start)

def plotProbability(vortex_prob, map, lower_bound, grid_spacing, tornado_track, title, file_name, obs=None, centers=None, min_prob=0.1, objects=None):
    nx, ny = vortex_prob.shape
    gs_x, gs_y = grid_spacing
    lb_x, lb_y = lower_bound
#   print lb_x, lb_y

    xs, ys = np.meshgrid(gs_x * np.arange(nx), gs_y * np.arange(ny))

    prob_color_map = matplotlib.cm.RdYlBu_r
    prob_color_map.set_under('#ffffff')

    track_xs, track_ys = map(*reversed(tornado_track))

    pylab.figure(figsize=(10, 8))
    pylab.axes((0.025, 0.025, 0.95, 0.925))

#   Plot Vortex Probability
    pylab.pcolormesh(xs, ys, vortex_prob, cmap=prob_color_map, vmin=min_prob, vmax=1.0)
    pylab.colorbar()#orientation='horizontal', aspect=40)

    pylab.plot(track_xs, track_ys, 'mv-', lw=2.5, mfc='k', ms=8)

    if centers:
        lines_x, lines_y = reorganizeCenters(centers, range(14400, 18300, 300))

        for line_x, line_y in zip(lines_x, lines_y):
            pylab.plot(line_x - lb_y, line_y - lb_x, 'ko-', markersize=2, linewidth=1)

    if obs is not None:
        pylab.contour(x, y, obs, levels=[0.95], colors='k')

    if objects:
        for obj in objects:
            obj_t, obj_y, obj_x = removeNones(obj, (13,) + vortex_prob.shape)

            if findSize(obj) >= 8:
                pylab.plot(1000 * np.array([obj_x.start, obj_x.start, obj_x.stop, obj_x.stop, obj_x.start]), 1000 * np.array([obj_y.start, obj_y.stop, obj_y.stop, obj_y.start, obj_y.start]), color='k', zorder=10)

    drawPolitical(map, scale_len=0) # scale_len=(xs[-1, -1] - xs[0, 0]) / 10000.)

    pylab.title(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def plotTiming(vortex_prob, vortex_times, times, map, grid_spacing, tornado_track, title, file_name, obs=None, centers=None, min_prob=0.1):
    nx, ny = vortex_prob.shape
    gs_x, gs_y = grid_spacing

    xs, ys = np.meshgrid(gs_x * np.arange(nx), gs_y * np.arange(ny))

    time_color_map = matplotlib.cm.Accent
    time_color_map.set_under('#ffffff')

    vortex_times = np.where(vortex_prob >= min_prob, vortex_times, -1)

    track_xs, track_ys = map(*reversed(tornado_track))

    pylab.figure(figsize=(10, 8))
    pylab.axes((0.025, 0.025, 0.95, 0.925))

    pylab.pcolormesh(xs, ys, vortex_times, cmap=time_color_map, vmin=times.min(), vmax=times.max())

    tick_labels = [ (datetime(2009, 6, 5, 18, 0, 0) + timedelta(seconds=int(t))).strftime("%H%M") for t in times ]
    bar = pylab.colorbar()#orientation='horizontal', aspect=40)
    bar.locator = FixedLocator(times)
    bar.formatter = FixedFormatter(tick_labels)
    bar.update_ticks()

    pylab.plot(track_xs, track_ys, 'mv-', lw=2.5, mfc='k', ms=8)

    drawPolitical(map, scale_len=(xs[-1, -1] - xs[0, 0]) / 10.)

    pylab.title(title)
    pylab.savefig(file_name)
    pylab.close()

def pickBox(tornado_track, objects):
    track_xs, track_ys = tornado_track

    track_xs = np.array(track_xs)
    track_ys = np.array(track_ys)

    keep_box = []
    for obj in objects:
        obj_t, obj_y, obj_x = removeNones(obj, (13, 255, 255))

        if np.all(track_xs >= obj_x.start * 1000) and np.all(track_xs < obj_x.stop * 1000) and np.all(track_ys >= obj_y.start * 1000) and np.all(track_ys < obj_y.stop * 1000):
            if keep_box == [] or findSize(obj) > find_size(keep_box):
                keep_box = obj

    return keep_box

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--parameter', dest='parameter', default='vort')
    ap.add_argument('--height', dest='interp_height', type=float, default=75.)
    ap.add_argument('--tag', dest='tag', required=True)
    ap.add_argument('--min-prob', dest='min_prob', type=float, default=0.1)
    ap.add_argument('--bbox', dest='bbox', action='store_true', default=False)

    args = ap.parse_args()

    exp_names = { 'prtrgn=1':"Storm Perturbations Only", 'newbc':"BC: $r_{0h}$ = 12 km", 'r0h=4km':'$r_{0h}$ = 4km', 'snd-no-w':'Sndgs do not update $w$', 'control':'Control ($r_{0h}$ = 6 km, Sndgs update $w$)', 
        'zupdtpt':'Control', 'z-no-05XP':"No MWR05XP Data", 'z-no-mm-05XP':"No MM or MWR05XP Data", 'z-no-mm':"No MM Data", 'z-no-v2':"No VORTEX2 data", 'z-no-snd':"No Sounding Data",
        'bc7dBZ,5ms':r"$r_{0h}$ = 4 km, BC: $\sigma_Z$ = 7 dBZ, $\sigma_{v_r}$ = 5 m s$^{-1}$", 'bcmult=1.03':"BC: Mult. Inflation Factor = 1.03",
        'r0h=6km-bc7dBZ,5ms':r"BC: $\sigma_Z$ = 7 dBZ, $\sigma_{v_r}$ = 5 m s$^{-1}$", '5dBZ,3ms-bc7dBZ,5ms':r"$\sigma_Z$ = 5 dBZ, $\sigma_{v_r}$ = 3 m s$^{-1}$", 'ebr':"Modified: Error, BC, and $r_{0h}$",
        'no-mm': 'No MM', 'mm':'MM', '05XP':'MM + MWR05XP', 'outer':"Outer Domain" }

    param = args.parameter
    interp_height = args.interp_height
    tag = args.tag

    if param == 'vort':
        description = "$\zeta$ > 0.0075 s$^{-1}$"
        domain_bounds = (slice(90, 160), slice(90, 160))

        max_refl_all, max_refl_da, max_refl_fcst = None, None, None

        try:
            vortex_centers = None
            vortex_centers = cPickle.load(open("vortex_centers_%s-%dm.pkl" % (tag, int(interp_height)), 'r'))
        except IOError:
            print "Can't find vortex center file for %s, %d m." % (tag, int(interp_height))
            vortex_centers = None

    elif param == 'refl':
        threshold = 40
        description = "$Z$ > %d dBZ" % threshold
        domain_bounds = (slice(None), slice(None))

#       max_refl_all, max_refl_da, max_refl_fcst = cPickle.load(open("max_obs_refl.pkl", 'r'))

#       max_refl_all  = np.where(max_refl_all  >= threshold, np.ones(max_refl_all.shape ), np.zeros(max_refl_all.shape ))
#       max_refl_da   = np.where(max_refl_da   >= threshold, np.ones(max_refl_da.shape  ), np.zeros(max_refl_da.shape  ))
#       max_refl_fcst = np.where(max_refl_fcst >= threshold, np.ones(max_refl_fcst.shape), np.zeros(max_refl_fcst.shape))

        max_refl_all, max_refl_da, max_refl_fcst = None, None, None

        vortex_centers = None

    elif param == 'w':
        threshold = 5.
        description = "$w$ > %d m s$^{-1}$" % threshold
        domain_bounds = (slice(None), slice(None))

        max_refl_all, max_refl_da, max_refl_fcst = None, None, None
        vortex_centers = None

    grid_spacing = goshen_1km_gs

    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs, domain_bounds[::-1])
    map = Basemap(**proj)

    max_prob = cPickle.load(open("max_%s_prob_%dm_%s.pkl" % (param, interp_height, tag),    'r'))
    argmax_prob = cPickle.load(open("argmax_%s_prob_%dm_%s.pkl" % (param, interp_height, tag), 'r'))

    tornado_track = zip(*((41.63,-104.383), (41.6134,-104.224)))

    if len(max_prob) >= 3:
        if param == 'w':
            max_prob_all, max_prob_da, max_prob_fcst, objects = max_prob
        else:
            max_prob_all, max_prob_da, max_prob_fcst = max_prob
            objects = None

        argmax_prob_all, argmax_prob_da, argmax_prob_fcst = argmax_prob
    else:
        if param == 'w':
            max_prob_all, objects = max_prob
        else:
            max_prob_all = max_prob[0]
            objects = None

        argmax_prob_all = argmax_prob[0]

    if param == 'w':
        bounding_box = pickBox(map(*reversed(tornado_track)), objects)
        cPickle.dump(bounding_box, open("bbox_%dm_%s.pkl" % (interp_height, tag), 'w'), -1)

    if args.bbox:
        bbox_buffer = 10
        bbox_offsets = [0, 10, 0] # t, y, x
        bbox = cPickle.load(open("bbox_2000m_%s.pkl" % tag, 'r'))
        objects = [ tuple(slice(b.start - bbox_buffer + o, b.stop + bbox_buffer + o) for b, o in zip(bbox, offsets)) ]

    start = datetime.utcnow()

    title = r"Probability of %s (%s at $z$ = %d m)" % (description, exp_names[tag], interp_height)
    lower_bounds = tuple([ b.start * g if b.start else 0 for b, g in zip(domain_bounds, grid_spacing) ])

    plotProbability(max_prob_all[domain_bounds],  map, lower_bounds, grid_spacing, tornado_track, 
        title, "max_%s_prob_%dm_%s_all.png" % (param, interp_height, tag), obs=max_refl_all, centers=vortex_centers, min_prob=args.min_prob, objects=objects)
    plotProbability(max_prob_da[domain_bounds],   map, lower_bounds, grid_spacing, tornado_track, 
        title, "max_%s_prob_%dm_%s_da.png" % (param, interp_height, tag), obs=max_refl_da, centers=vortex_centers, min_prob=args.min_prob)
    plotProbability(max_prob_fcst[domain_bounds], map, lower_bounds, grid_spacing, tornado_track, 
        title, "max_%s_prob_%dm_%s_fcst.png" % (param, interp_height, tag), obs=max_refl_fcst, centers=vortex_centers, min_prob=args.min_prob)

    title = r"Earliest Timing of %s (%s at $z$ = %d m)" % (description, exp_names[tag], interp_height)

    plotTiming(max_prob_all[domain_bounds],  argmax_prob_all[domain_bounds], np.arange(10800, 18300, 300),  map, grid_spacing, tornado_track, 
        title, "argmax_%s_prob_%dm_%s_all.png" % (param, interp_height, tag), obs=max_refl_all, centers=vortex_centers, min_prob=args.min_prob)
    plotTiming(max_prob_da[domain_bounds],   argmax_prob_da[domain_bounds], np.arange(10800, 14700, 300),   map, grid_spacing, tornado_track, 
        title, "argmax_%s_prob_%dm_%s_da.png" % (param, interp_height, tag), obs=max_refl_da, centers=vortex_centers, min_prob=args.min_prob)
    plotTiming(max_prob_fcst[domain_bounds], argmax_prob_fcst[domain_bounds], np.arange(14400, 18300, 300), map, grid_spacing, tornado_track, 
        title, "argmax_%s_prob_%dm_%s_fcst.png" % (param, interp_height, tag), obs=max_refl_fcst, centers=vortex_centers, min_prob=args.min_prob)

    print "Time to plot images:", datetime.utcnow() - start
    return

if __name__ == "__main__":
    main()
