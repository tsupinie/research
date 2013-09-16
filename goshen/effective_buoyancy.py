
import numpy as np
import scipy.spatial.distance as distance

import matplotlib
matplotlib.use('agg')
import pylab

import Nio as nio

from dataload import loadEnsemble
from temporal import goshen_1km_temporal
from grid import goshen_1km_grid
from computeQuantities import theta2Temperature
from util import decompressVariable

def computeDensity(**kwargs):
    t = theta2Temperature(**kwargs)
    return kwargs['p'] / (t * 287.)

def horiz_laplacian(field, coords):
    ys, xs = coords
    d2field_dx2 = np.empty(field.shape)
    d2field_dy2 = np.empty(field.shape)

    d2field_dx2[:, :, 1:-1] = (field[:, :, 2:] - 2 * field[:, :, 1:-1] + field[:, :, :-2]) / (xs[:, :, 1:-1] - xs[:, :, :-2]) ** 2
    d2field_dx2[:, :, 0]    = d2field_dx2[:, :, 1]
    d2field_dx2[:, :, -1]   = d2field_dx2[:, :, -2]

    d2field_dy2[:, 1:-1, :] = (field[:, 2:, :] - 2 * field[:, 1:-1, :] + field[:, :-2, :]) / (ys[:, 1:-1, :] - ys[:, :-2, :]) ** 2
    d2field_dy2[:, 0, :]    = d2field_dy2[:, 1, :]
    d2field_dy2[:, -1, :]   = d2field_dy2[:, -2, :]

    return d2field_dx2 + d2field_dy2

def computeCellVolume(coords):
    zs, ys, xs = coords
    volume = np.empty(zs.shape)

    zs_int = (zs[:-1, :, :] + zs[1:, :, :]) / 2
    ys_int = (ys[:, :-1, :] + ys[:, 1:, :]) / 2
    xs_int = (xs[:, :, :-1] + xs[:, :, 1:]) / 2

    volume[1:-1, 1:-1, 1:-1] = (zs_int[1:, 1:-1, 1:-1] - zs_int[:-1, 1:-1, 1:-1]) * (ys_int[1:-1, 1:, 1:-1] - ys_int[1:-1, :-1, 1:-1]) * (xs_int[1:-1, 1:-1, 1:] - xs_int[1:-1, 1:-1, :-1])
    volume[0, :, :], volume[-1, :, :] = volume[1, :, :], volume[-2, :, :]
    volume[:, 0, :], volume[:, -1, :] = volume[:, 1, :], volume[:, -2, :]
    volume[:, :, 0], volume[:, :, -1] = volume[:, :, 1], volume[:, :, -2]
    return volume

def computeSurfaceNormals(coords):
    return

def effectiveBuoyancy(density, coords, plane=None):
    g = 9.806

    laplacian_density = horiz_laplacian(density, coords[1:])

    zs, ys, xs = coords
    volume = computeCellVolume(coords)

#   coord_array = np.vstack((zs.ravel(), ys.ravel(), xs.ravel())).T

    cur_z = -1
    from datetime import timedelta, datetime

    dist_time = timedelta(0)
    buoy_time = timedelta(0)

    total_start = datetime.now()
    box_size = 15

    if plane is None:
        comp_bounds = (slice(None), slice(None), slice(None))
        lower_bound = (0, 0, 0)
    else:
        comp_bounds = ()
        lower_bound = ()
        for ax in [ 'z', 'y', 'x' ]:
            if ax in plane:
                comp_bounds = comp_bounds + (plane[ax],)
                lower_bound = lower_bound + (plane[ax],)
            else:
                comp_bounds = comp_bounds + (slice(None),)
                lower_bound = lower_bound + (0,)

    eff_buoy = np.empty(np.transpose(np.atleast_3d(density[comp_bounds]), (2, 0, 1)).shape)

    try:
        for idx in np.ndindex(eff_buoy.shape):
            if idx[0] != cur_z:
                print "Working on level %d" % (idx[0] + 1)
                cur_z = idx[0]

            lb_z, lb_y, lb_x = lower_bound

            z, y, x = zs[idx], ys[idx], xs[idx]
            bound_z = slice(max(0, idx[0] + lb_z - box_size), min(zs.shape[0], idx[0] + lb_z + box_size))
            bound_y = slice(max(0, idx[1] + lb_y - box_size), min(ys.shape[1], idx[1] + lb_y + box_size))
            bound_x = slice(max(0, idx[2] + lb_x - box_size), min(xs.shape[2], idx[2] + lb_x + box_size))

            bounds = (bound_z, bound_y, bound_x)

#           point = np.array([zs[idx], ys[idx], xs[idx]])[np.newaxis, :]
            start = datetime.now()
            dists = np.sqrt((z - zs[bounds]) ** 2 + (y - ys[bounds]) ** 2 + (x - xs[bounds]) ** 2)
#           dists = distance.cdist(coord_array, point).reshape(density.shape)
            dist_time += (datetime.now() - start)

            start = datetime.now()
            integrand = (laplacian_density[bounds] * volume[bounds] / dists)
            integrated = integrand[np.where(np.isfinite(integrand))].sum()
            eff_buoy[idx] = g / (4 * np.pi) * integrated
            buoy_time += (datetime.now() - start)

    except KeyboardInterrupt:
        print "Total time:", datetime.now() - total_start
        print "Time to compute distances:", dist_time
        print "Time to compute buoyancy:", buoy_time
        import sys
        sys.exit()

    print "Total time:", datetime.now() - total_start
    print "Time to compute distances:", dist_time
    print "Time to compute buoyancy:", buoy_time

    return eff_buoy

def main():
    base_path = "/caps2/tsupinie/1kmf-control/"
    temp = goshen_1km_temporal(start=14400, end=14400)
    grid = goshen_1km_grid()
    n_ens_members = 40

    np.seterr(all='ignore')

    ens = loadEnsemble(base_path, [ 11 ], temp.getTimes(), ([ 'pt', 'p' ], computeDensity))
    ens = ens[0, 0]

    zs = decompressVariable(nio.open_file("%s/ena001.hdfgrdbas" % base_path, mode='r', format='hdf').variables['zp'])
    xs, ys = grid.getXY()
    xs = xs[np.newaxis, ...].repeat(zs.shape[0], axis=0)
    ys = ys[np.newaxis, ...].repeat(zs.shape[0], axis=0)

    eff_buoy = effectiveBuoyancy(ens, (zs, ys, xs), plane={'z':10})
    print eff_buoy

    pylab.figure()
    pylab.contourf(xs[0], ys[0], eff_buoy[0], cmap=matplotlib.cm.get_cmap('RdBu_r'))
    pylab.colorbar()

    grid.drawPolitical()

    pylab.suptitle("Effective Buoyancy")
    pylab.savefig("eff_buoy.png")
    pylab.close()
    return

if __name__ == "__main__":
    main()
