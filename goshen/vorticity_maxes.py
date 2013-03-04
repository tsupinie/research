
import numpy as np
import scipy.ndimage as ndimage

import Nio as nio

import glob
import cPickle

from util import loadAndInterpolateEnsemble, setupMapProjection
from computeQuantities import computeVorticity, computeReflectivity

import pylab
from mpl_toolkits.basemap import Basemap

def plotCenters(vortex_centers, proj, title, file_name, n_ens_target=None, t_ens_target=None, vort=None, vectors=None):
    pylab.figure()
    map = Basemap(**proj)

    if vort is not None:
        xs, ys = np.meshgrid(1000 * np.arange(vort.shape[1]), 1000 * np.arange(vort.shape[0]))
        pylab.contourf(xs, ys, vort, levels=np.arange(-0.02, 0.025, 0.005))
        pylab.colorbar()

    if vectors is not None:
        xs, ys = np.meshgrid(1000 * np.arange(vort.shape[1]), 1000 * np.arange(vort.shape[0]))
        stride = 2
        bounds = (slice(None, None, stride), slice(None, None, stride))

        vector_shape = vectors.shape
        us, vs = zip(*vectors.ravel())
        us = np.array(us).reshape(vector_shape)
        vs = np.array(vs).reshape(vector_shape)

        pylab.quiver(xs[bounds], ys[bounds], us[bounds], vs[bounds])

    for (n_ens, t_ens), center in vortex_centers.iteritems():
        if (n_ens_target is not None and n_ens_target != n_ens) or (t_ens_target is not None and t_ens_target != t_ens):
            continue

        x, y = center
        map.plot(x, y, 'ko')
#       pylab.text(x, y, n_ens, va='bottom', ha='right')

    map.readshapefile("countyp020", 'counties', linewidth=0.5)
    map.drawstates(linewidth=1.0)
    map.drawcoastlines(linewidth=1.0)
    map.drawcountries(linewidth=1.0)

    pylab.title(title)
    pylab.savefig(file_name)
    return

def main():
    interp_height = 975
    vort_thresh = 0.01
    var_list = ['u', 'v', 'dx', 'dy']

    grdbas = nio.open_file('1km-control/ena001.hdfgrdbas', mode='r', format='hdf')
    xs = grdbas.variables['x'][:]
    ys = grdbas.variables['y'][:]

    x_grid_spacing = xs[1] - xs[0]
    y_grid_spacing = ys[1] - ys[0]
    nx = xs.shape[0]
    ny = ys.shape[0]

    files = glob.glob("/caps1/tsupinie/1km-control-20120712/ena???.hdf014700")

#   bounds = (slice(0, nx), slice(0, ny))
    bounds = (slice(90, 160), slice(90, 160))
    proj = {'projection':'lcc', 'resolution':'l',
        'width':(x_grid_spacing * nx), 'height':(y_grid_spacing * ny),
        'ctr_lat':41.61975, 'ctr_lon':-104.34843, 
        'lat_0':40.61998, 'lon_0':-107.344, 'lat_1':30., 'lon_1':60.
    }
    proj = setupMapProjection(proj, (x_grid_spacing, y_grid_spacing), bounds=bounds)

    def vectorTuple(**kwargs):
        vectors = np.empty(kwargs['u'].shape, dtype=[('u', np.float32), ('v', np.float32)])
        vectors['u'] = kwargs['u']
        vectors['v'] = kwargs['v']
        return vectors

    ens_vort, ens_members, ens_times = loadAndInterpolateEnsemble(files, var_list, computeVorticity, "1km-control/ena001.hdfgrdbas", points={'z':interp_height})
    ens_vectors, ens_members, ens_times = loadAndInterpolateEnsemble(files, ['u', 'v'], vectorTuple, "1km-control/ena001.hdfgrdbas", points={'z':interp_height})

    ens_vort = ens_vort[:, :, :-2, :-2]
    ens_vectors = ens_vectors[:, :, :-2, :-2]

    correct_region = set()
    correct_members = [ "%03d" % m for m in [1, 4, 5, 6, 7, 11, 12, 14, 17, 19, 21, 23, 24, 25, 26, 27, 30, 31, 33, 35, 36, 38, 40] ]

    vort_centers = {}
    for lde in range(len(ens_members)):
        for wdt in range(len(ens_times)):
            this_vort = ens_vort[lde, wdt]

            vort_mask = np.where(this_vort > vort_thresh, True, False)
            vort_maxes, num_vort_maxes = ndimage.label(vort_mask)
            objects = ndimage.measurements.find_objects(vort_maxes)

#           print objects

            if len(objects) > 0:
                largest_vort_max_idx =  np.argmax([len(this_vort[o].ravel()) for o in objects])
            else:
                print "n_ens = %s, t_ens = %06d has no areas of vorticity > %f" % (ens_members[lde], ens_times[wdt], vort_thresh)
                continue

            region = objects[largest_vort_max_idx]

            if ens_members[lde] in correct_members:
                y_idxs, x_idxs = np.meshgrid(*[ np.arange(c.start, c.stop, dtype=int) for c in region ])
                correct_region = correct_region | set(zip(y_idxs.ravel(), x_idxs.ravel()))

            largest_vort_max = this_vort[region]
            center = zip(*np.where(this_vort == largest_vort_max.max()))[0]
            center = tuple(reversed(center))
            vort_centers[(ens_members[lde], ens_times[wdt])] = tuple([ (c - b.start) * gs for c, gs, b in zip(center, [x_grid_spacing, y_grid_spacing], bounds) ])

    correct_lb_y, correct_lb_x = [ min(p) for p in zip(*correct_region) ] 
    correct_ub_y, correct_ub_x = [ max(p) for p in zip(*correct_region) ] 

    correct_members = [ "%03d" % m for m in [1, 2, 4, 5, 6, 7, 9, 11, 12, 14, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 35, 36, 37, 38, 40] ]
    for lde in range(len(ens_members)):
        for wdt in range(len(ens_times)):
            this_vort = ens_vort[lde, wdt]
            object = this_vort[correct_lb_y:correct_ub_y, correct_lb_x:correct_ub_x]

            center = zip(*np.where(this_vort == object.max()))[0]
            center = tuple(reversed(center))

            x, y = center
            if ens_members[lde] in [ '003', '010', '013', '015', '020', '039' ]:
                object = this_vort[(y + 2):(y + 10), correct_lb_x:correct_ub_x]
            elif ens_members[lde] in [ '008', '016' ]:
                object = this_vort[correct_lb_y:correct_ub_y, (x - 10):(x - 2)]
            elif ens_members[lde] in [ '034' ]:
                object = this_vort[(y - 10):(y - 2), correct_lb_x:correct_ub_x]
            elif ens_members[lde] in [ '032' ]:
                del vort_centers[(ens_members[lde], ens_times[wdt])]

            print ens_members[lde], ens_times[wdt], object.max(), object

            if np.all(object < vort_thresh) and (ens_members[lde], ens_times[wdt]) in vort_centers:
                del vort_centers[(ens_members[lde], ens_times[wdt])]
                print "n_ens = %s, t_ens = %06d relocated vorticity max < %f, deleting" % (ens_members[lde], ens_times[wdt], vort_thresh)
                continue

            center = zip(*np.where(this_vort == object.max()))[0]
            center = tuple(reversed(center))

            vort_centers[(ens_members[lde], ens_times[wdt])] = tuple([ (c - b.start) * gs for c, gs, b in zip(center, [x_grid_spacing, y_grid_spacing], bounds) ])

    cPickle.dump(vort_centers, open('vort-centers-20120712.pkl', 'w'), -1)

    plotCenters(vort_centers, proj, "Vortex Centers at 05 June 2009 2205 UTC", "vort_centers.png", t_ens_target=14700, vort=ens_vort.max(axis=0).max(axis=0)[tuple(reversed(bounds))])

    for lde, n_ens in enumerate(ens_members):
        this_vort = ens_vort[lde, 0]
        this_vector = ens_vectors[lde, 0]
        plotCenters(vort_centers, proj, "Vortex Center (Member %s) at 05 June 2009 2205 UTC" % n_ens, "vort_center_%s.png" % n_ens, n_ens_target=n_ens, t_ens_target=14700, vort=this_vort[tuple(reversed(bounds))], vectors=this_vector[tuple(reversed(bounds))])

if __name__ == "__main__":
    main()
