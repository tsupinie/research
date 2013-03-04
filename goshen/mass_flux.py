
from util import loadAndInterpolateEnsemble
from computeQuantities import theta2Temperature

import matplotlib
matplotlib.use('agg')
import pylab
import numpy as np

import glob, argparse, cPickle

def getWAndDensity(**kwargs):
    wrho = np.empty(kwargs['w'].shape, dtype=[('w', np.float32), ('rho', np.float32)])
    temp = theta2Temperature(**kwargs)

    wrho['w'] = kwargs['w']
    wrho['rho'] = kwargs['p'] / (temp * 287.)
    return wrho

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp-name', dest='exp_name', required=True)

    args = ap.parse_args()

    base_path = "/caps1/tsupinie/1km-control-%s" % args.exp_name

    files = glob.glob("%s/ena???.hdf014[47]00" % base_path)
    files.extend(glob.glob("%s/ena???.hdf01[5678]?00" % base_path))

    ens_wrho, ens_members, ens_times = loadAndInterpolateEnsemble(files, ['w', 'pt', 'p'], getWAndDensity, "%s/ena001.hdfgrdbas" % base_path)

    print ens_wrho.shape
    w_mean = np.maximum(0, ens_wrho['w'].mean(axis=0))
    rho_mean = ens_wrho['rho'].mean(axis=0)

    flux_box = None
    box_size = 10
    w_test_level = 12

    track_cutoff = 16200
    dt_ens = ens_times[1] - ens_times[0]
    num_track_steps = 1 + (track_cutoff - ens_times[0]) / dt_ens
   
    flux_shape = list(w_mean.shape[:2])
    flux_shape.extend([ 2 * box_size + 1 ] * 2)
    mass_flux = np.empty(flux_shape, dtype=w_mean.dtype)

    first_flux_box = None
    last_flux_box = None

    for wdt, t_ens in enumerate(ens_times):
        if flux_box is None:
            # Create the first box, keep track of this first position
            w_test_data = w_mean[wdt, w_test_level]
            w_max_y, w_max_x = np.unravel_index(np.argmax(w_test_data), w_test_data.shape)
            flux_box = (slice(None), slice(w_max_y - box_size, w_max_y + box_size + 1), slice(w_max_x - box_size, w_max_x + box_size + 1))
            first_flux_box = flux_box
        elif t_ens <= track_cutoff:
            # Track the updraft through the cutoff time, keep track of the last position of the box
            w_test_data = w_mean[wdt, w_test_level][flux_box[1:]]
            old_box_lb_y, old_box_lb_x = [ b.start for b in flux_box[1:] ]
            w_max_y, w_max_x = np.unravel_index(np.argmax(w_test_data), w_test_data.shape)

            w_max_y += old_box_lb_y
            w_max_x += old_box_lb_x

            flux_box = (slice(None), slice(w_max_y - box_size, w_max_y + box_size + 1), slice(w_max_x - box_size, w_max_x + box_size + 1))
            last_flux_box = flux_box
        else:
            # Extrapolate the box for the remainder of the forecast, with a fudge factor
            fudge_factor = 4

            steps_since_cutoff = (t_ens - track_cutoff) / dt_ens

            first_lb_y, first_lb_x = [ b.start for b in first_flux_box[1:] ]
            last_lb_y, last_lb_x = [ b.start for b in last_flux_box[1:] ]
            
            per_step_trans_x = (last_lb_x - first_lb_x) / num_track_steps
            per_step_trans_y = (last_lb_y - first_lb_y) / num_track_steps

            box_lb_x = per_step_trans_x * steps_since_cutoff + last_lb_x
            box_lb_y = (per_step_trans_y + fudge_factor) * steps_since_cutoff + last_lb_y

            flux_box = (slice(None), slice(box_lb_y, box_lb_y + 2 * box_size + 1), slice(box_lb_x, box_lb_x + 2 * box_size + 1))

        print flux_box[1:]

        mass_flux[wdt] = w_mean[wdt][flux_box] * (1000 * 1000) * rho_mean[wdt][flux_box]

    total_mass_flux = mass_flux.sum(axis=-1).sum(axis=-1)
    print total_mass_flux.shape

    cPickle.dump(total_mass_flux, open("mass_flux_%s.pkl" % args.exp_name, 'w'), -1)
    return

if __name__ == "__main__":
    main()
