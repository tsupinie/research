
import Nio as nio

import numpy as np
from scipy.interpolate import interp1d
from scipy import weave
import matplotlib
matplotlib.use('agg')
import pylab

import argparse
import glob
import cPickle
import gc
from datetime import datetime, timedelta

from computeQuantities import computeVorticity, computeReflectivity
from dataload import loadEnsemble

def interpolate_old(data, point, axes):
    interp_grid = data
    for ax in range(len(data.shape) - 1, -1, -1):
        if ax in point:
            if len(axes[ax].shape) == 1:
                interp_grid = interp1d(axes[ax], interp_grid, axis=ax, kind='linear')(point[ax])
            else:
                interp_ax = axes[ax]
                for ax_ax in range(len(data.shape) - 1, ax, -1):
                    if ax_ax in point:
                        interp_ax = interp1d(axes[ax_ax], interp_ax, axis=ax_ax, kind='linear')(point[ax_ax])

                new_dims = tuple([ d for nd, d in enumerate(interp_grid.shape) if nd != ax])
                temp_grid = np.empty(new_dims)
                for index in np.ndindex(new_dims):
                    slice_dims = list(index)
                    slice_dims.insert(ax, slice(None))
                    slice_dims = tuple(slice_dims)

                    temp_grid[index] = interp1d(interp_ax[slice_dims], interp_grid[slice_dims], axis=ax, kind='linear', bounds_error=False)(point[ax])
                interp_grid = temp_grid
    return interp_grid

def thresholdField(field, threshold_value):
    return np.where(field >= threshold_value, np.ones(field.shape), np.zeros(field.shape))

def generateFakeEnsemble(grid_dims, grid_spacing, vortex_center, nens, spread):
    def velocity(radius):
        max_v_theta = 50
        r_max_v_theta = 2 * grid_spacing

        if radius <= r_max_v_theta:
            return max_v_theta * radius / r_max_v_theta
        else:
            return max_v_theta * r_max_v_theta / radius

    nx, ny = grid_dims

    u = np.empty((nens, ny, nx))
    v = np.empty((nens, ny, nx))

    vortex_center_x, vortex_center_y = vortex_center

    vortex_center_xs = np.random.normal(vortex_center_x, spread, nens)
    vortex_center_ys = np.random.normal(vortex_center_y, spread, nens)

    for nen, jdy, idx in np.ndindex(*u.shape):
        delta_x = grid_spacing * idx - vortex_center_xs[nen]
        delta_y = grid_spacing * jdy - vortex_center_ys[nen]
        theta = np.arctan2(delta_y, delta_x)
        tan_vel = velocity(np.hypot(delta_x, delta_y))

        u[nen, jdy, idx] = -tan_vel * np.sin(theta)
        v[nen, jdy, idx] = tan_vel * np.cos(theta)
    return u, v   

def r_enumerate(L):
    for idx in xrange(len(L) - 1, -1, -1):
        yield idx, L[idx]

def swath(data, threshold, nens, times, lower_p_bound):
    time_integrated_data = data.max(axis=1)
    ti_prob = thresholdField(time_integrated_data, threshold).sum(axis=0) / nens

    prob = thresholdField(data, threshold).sum(axis=0) / nens
    max_prob = prob.max(axis=0)
    argmax_prob = np.argmax(prob, axis=0)

    for wdt, time in r_enumerate(times):
        argmax_prob = np.where(argmax_prob == wdt, time, argmax_prob)

#   argmax_prob = np.where(max_prob < lower_p_bound, -1, argmax_prob)
    return ti_prob, argmax_prob

def findProbObjects(data, data_threshold, prob_threshold):
    from scipy.ndimage.measurements import label, find_objects
    connectivity = np.ones((3, 3, 3), dtype=np.int32)

    ens_prob = (data >= data_threshold).sum(axis=0) / float(data.shape[0])
    labels, n_objs = label(ens_prob >= prob_threshold, structure=connectivity)
    return find_objects(labels)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-dir', dest='data_dir', default=None)
    ap.add_argument('--parameter', dest='parameter', default='vort')
    ap.add_argument('--height', dest='interp_height', type=float, default=75.)
    ap.add_argument('--tag', dest='tag', required=True)

    args = ap.parse_args()

    if args.data_dir is None:
        nx, ny = 50, 50
        grid_spacing = 1000
        n_ensemble_members = 40
        times = range(0, 1800, 300)

        u_all = np.empty((len(times), n_ensemble_members, ny, nx))
        v_all = np.empty((len(times), n_ensemble_members, ny, nx))

        for wdt, t_ens in enumerate(times):
            u, v = generateFakeEnsemble((nx, ny), grid_spacing, (grid_spacing * nx / 2, grid_spacing * (ny / 2 - (900 - t_ens) / 90.)), n_ensemble_members, (2 + t_ens / 900.) * grid_spacing)

            u_all[wdt] = u
            v_all[wdt] = v

        u_all = np.transpose(u_all, (1, 0, 2, 3))
        v_all = np.transpose(v_all, (1, 0, 2, 3))
    else:
#       files = glob.glob("%s/ena???.hdf0*" % args.data_dir)
#       files = glob.glob("%s/ena???.hdf010800" % args.data_dir)

        if args.parameter == 'vort':
            var_list = ['u', 'v', 'dx', 'dy']
            func = computeVorticity
        elif args.parameter == 'refl':
            func = computeReflectivity
            var_list = ['p', 'pt', 'qr', 'qs', 'qh']
        elif args.parameter == 'w':
            func = lambda **kwargs: kwargs['w']
            var_list = [ 'w' ]

        n_ensemble_members = 40
        times = np.arange(14400, 18300, 300)
        param_all = loadEnsemble(args.data_dir, n_ensemble_members, times, (var_list, func), { 'z':args.interp_height }, agl=True)

    cutoff = np.where(times == 14400)[0]

    if args.parameter == 'vort':
        threshold = 0.0075
        lower_p_bound = 0.2
    elif args.parameter == 'refl':
        threshold = 40.
        lower_p_bound = 0.1
    elif args.parameter == 'w':
        threshold = 5.
        lower_p_bound = 0.1

        objects = findProbObjects(param_all, threshold, 0.2)

    max_prob_all, argmax_prob_all = swath(param_all, threshold, n_ensemble_members, times, lower_p_bound)

    if cutoff < len(times):
        max_prob_da,   argmax_prob_da   = swath(param_all[:, :(cutoff + 1), :, :], threshold, n_ensemble_members, times[:(cutoff + 1)], lower_p_bound)
        max_prob_fcst, argmax_prob_fcst = swath(param_all[:, cutoff:, :, :],       threshold, n_ensemble_members, times[cutoff:],       lower_p_bound)
 
        if args.parameter == 'w':
            cPickle.dump((max_prob_all,    max_prob_da,    max_prob_fcst,   objects), open("max_%s_prob_%dm_%s.pkl" % (args.parameter, args.interp_height, args.tag),    'w'), -1)
            cPickle.dump((argmax_prob_all, argmax_prob_da, argmax_prob_fcst        ), open("argmax_%s_prob_%dm_%s.pkl" % (args.parameter, args.interp_height, args.tag), 'w'), -1)
        else:
            cPickle.dump((max_prob_all,    max_prob_da,    max_prob_fcst),    open("max_%s_prob_%dm_%s.pkl" % (args.parameter, args.interp_height, args.tag),    'w'), -1)
            cPickle.dump((argmax_prob_all, argmax_prob_da, argmax_prob_fcst), open("argmax_%s_prob_%dm_%s.pkl" % (args.parameter, args.interp_height, args.tag), 'w'), -1)
    else:
        if args.parameter == 'w':
            cPickle.dump((max_prob_all, objects),    open("max_%s_prob_%dm_%s.pkl" % (args.parameter, args.interp_height, args.tag),    'w'), -1)
            cPickle.dump((argmax_prob_all,),         open("argmax_%s_prob_%dm_%s.pkl" % (args.parameter, args.interp_height, args.tag), 'w'), -1)
        else:
            cPickle.dump((max_prob_all,),    open("max_%s_prob_%dm_%s.pkl" % (args.parameter, args.interp_height, args.tag),    'w'), -1)
            cPickle.dump((argmax_prob_all,), open("argmax_%s_prob_%dm_%s.pkl" % (args.parameter, args.interp_height, args.tag), 'w'), -1)

    return

if __name__ == "__main__":
    main()
