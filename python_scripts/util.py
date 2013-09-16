"""
util.py
Author:         Tim Supinie (tsupinie@ou.edu)
Last Modified:  21 February 2013
Purpose:        A bunch of helper functions for my research.  Stuff like loading and interpolating data and plotting maps.
"""

import numpy as np

from matplotlib import transforms

import copy
from datetime import datetime, timedelta
import cPickle
from math import ceil
from multiprocessing import Process, Queue
import string

_interp_code = ""
_lib_code = ""


def decompressVariable(hd_var, dindex=None):
    if not hasattr(hd_var, 'min') or not hasattr(hd_var, 'max'):
        if dindex is None:
            return hd_var[:]
        else:
            return hd_var[dindex]

    max_i16 = np.iinfo(np.int16).max
    min_i16 = np.iinfo(np.int16).min

    if dindex is None:
        fraction = (hd_var[:].astype(np.float32) - min_i16) / (max_i16 - min_i16)

        new_shape = [ hd_var.shape[0] ]
        for idx in range(len(hd_var.shape) - 1): new_shape.append(1)
        new_shape = tuple(new_shape)

        return hd_var.min.reshape(new_shape) * (1 - fraction) + hd_var.max.reshape(new_shape) * (fraction)
    else:
        fraction = (hd_var[dindex].astype(np.float32) - min_i16) / (max_i16 - min_i16)

        if type(dindex[0]) == int:
            new_shape = [ 1 ]
        else:
            start = dindex[0].start
            if dindex[0].start is None: start = 0

            stop = dindex[0].stop 
            if dindex[0].stop is None: stop = hd_var.shape[0]

            step = dindex[0].step 
            if dindex[0].step is None: step = 1

            new_shape = [ len(range(start, stop, step)) ]

        for idx in range(len(fraction.shape) - 1): new_shape.append(1)
        new_shape = tuple(new_shape)

        decompressed_data = hd_var.min[dindex[0]].reshape(new_shape) * (1 - fraction) + hd_var.max[dindex[0]].reshape(new_shape) * (fraction)
        return decompressed_data

def isolateObsTimes(obs, times, round=True):
    epoch = datetime(1970, 1, 1, 0, 0, 0)

    dtype = sorted(obs.dtype.fields.items(), key=lambda (f, (d, o)): o)
    new_dtype = [ (f, d) for (f, (d, o)) in dtype ]
    new_dtype.append(('nom_time', np.dtype('float64')))

    keep_obs = [ ] #np.empty((0,), dtype=new_dtype)

    too_early = 0
    too_late = 0

    for idx in xrange(obs.shape[0]):
        ob = obs[idx]
        ob_time = epoch + timedelta(seconds=ob['time'])

        if ob_time in times:
            keep_obs.append(tuple(ob) + (ob['time'],))
        elif ob_time >= min(times) and ob_time <= max(times): # and ob['latitude'] not in keep_obs['latitude'] and ob['longitude'] not in keep_obs['longitude']:
            if round:
                round_time = np.round(ob['time'] / 300) * 300
                keep_obs.append(tuple(ob) + (round_time,))

        elif ob_time < min(times):
            too_early += 1
        elif ob_time > max(times):
            too_late += 1

#   print "# obs too early:", too_early
#   print "# obs too late:", too_late
#   print keep_obs.shape

    return np.array(keep_obs, dtype=new_dtype)

def thinObs(obs, map, width, height):
    horiz_distance_thresh = 1000
    pressure_thresh = 5
    time_thresh = 300

    obs_xs, obs_ys = map(obs['longitude'], obs['latitude'])

    start_idx = 0
    while obs_xs[start_idx] < 0 or obs_xs[start_idx] > width or obs_ys[start_idx] < 0 or obs_ys[start_idx] > height:
        start_idx += 1

    keep_obs = np.array([ obs[start_idx] ], dtype=obs.dtype)
    keep_xs = np.array([ obs_xs[start_idx] ])
    keep_ys = np.array([ obs_ys[start_idx] ])
    for ob, x, y in zip(obs[(start_idx + 1):], obs_xs[(start_idx + 1):], obs_ys[(start_idx + 1):]):
        if x >= 0 and x <= width and y >= 0 and y <= height:

            horiz_dists = np.hypot(keep_xs - x, keep_ys - y)
            pres_dists = np.abs(keep_obs['pres'] - ob['pres'])
            time_dists = np.abs(keep_obs['nom_time'] - ob['nom_time'])

            if np.where((horiz_dists < horiz_distance_thresh) & (pres_dists < pressure_thresh) & (time_dists < time_thresh))[0].shape[0] == 0:
                keep_obs.resize(keep_obs.shape[0] + 1)
                keep_obs[-1] = tuple(ob)

                keep_xs.resize(keep_xs.shape[0] + 1)
                keep_xs[-1] = x

                keep_ys.resize(keep_ys.shape[0] + 1)
                keep_ys[-1] = y

    return np.array(keep_obs, dtype=obs.dtype)

def rearrangeSoundingObs(sounding_obs):
    dtype = [('id', np.dtype('|S4')), ('obtype', np.dtype('|S8')), ('time', np.dtype('float64')), ('latitude', np.dtype('float64')), ('longitude', np.dtype('float64')), ('elevation', np.dtype('float64')), 
                ('temp', np.dtype('float64')), ('dewp', np.dtype('float64')), ('pres', np.dtype('float64')), ('wind_dir', np.dtype('float64')), ('wind_spd', np.dtype('float64'))]

    rearr_sounding_obs = []
    for sounding in sounding_obs:
        release_time = datetime.strptime(sounding['release_time'], "%Y, %m, %d, %H:%M:%S")
        release_epoch = (release_time - datetime(1970, 1, 1, 0, 0, 0)).total_seconds()
        for time, lat, lon, hght, temp, dewp, pres, wdir, wspd in zip(*[ sounding[k] for k in ['time', 'latitude', 'longitude', 'altitude', 'temperature', 'dewpoint', 'pressure', 'wind_dir', 'wind_spd']]):
            if temp != 999.0 and dewp != 999.0 and pres != 9999.0 and wdir != 999.0 and wspd != 999.0:
                temp_F = 9. / 5. * temp + 32 
                dewp_F = 9. / 5. * dewp + 32
                rearr_sounding_obs.append((sounding['release_site'][:4], "SNDG", ceil(release_epoch / 3600.) * 3600, lat, lon, hght, temp_F, dewp_F, pres, wdir, wspd))

    return np.array(rearr_sounding_obs, dtype=dtype)

def loadObs(file_names, times, map, dimensions, sounding_obs=None, round_time=True):
    """
    loadObs()
    Purpose:    Load all the observation data, throw out those observations outside the grid, thin the observation data, and put everything
                    in a single numpy record array.  The exact mechanics of this function probably don't really apply because your stuff will
                    be different, but it should return the same general thing.
    Parameters: [ Perhaps not important here.  Ask if you're interested. ]
    Returns:    A numpy record array where each entry in the record array is an observation.  The record array should have the following dtype:

                  [('id', np.dtype('|S4')), ('obtype', np.dtype('|S8')), ('time', np.dtype('float64')), ('latitude', np.dtype('float64')), 
                    ('longitude', np.dtype('float64')), ('elevation', np.dtype('float64')), ('temp', np.dtype('float64')), 
                    ('dewp', np.dtype('float64')), ('pres', np.dtype('float64')), ('wind_dir', np.dtype('float64')), 
                    ('wind_spd', np.dtype('float64'))]
    """
    obs = []

    for file_name in file_names:
        file_obs = cPickle.load(open(file_name, 'r'))

        if sounding_obs is not None and file_name in sounding_obs:
            file_obs = rearrangeSoundingObs(file_obs)

        obs.append(file_obs)

    all_obs = np.concatenate(tuple(obs)).astype(obs[0].dtype)

    keep_obs = isolateObsTimes(all_obs, times, round=round_time)
    keep_obs.sort(order='nom_time')

    keep_obs = thinObs(keep_obs, map, *dimensions)

    print keep_obs.shape
#   print keep_obs

    return keep_obs

def runConcurrently(target, placeholder_vals, args=[], kwargs={}, max_concurrent=-1, zip_result=False):
    """
    runConcurrently()
    Purpose:    Runs several instances of a function at the same time and returns all their outputs as a list.
    Parameters: target [type=function]
                    Function to run.  For this version of the function, it must return something.
                placeholder_vals [type=list,tuple]
                    Values of a placeholder parameter to run the function on.
                args [optional, type=list]
                    Arguments to pass to the target function.  One or more may have the special value "__placeholder__", which
                    will be replaced with a value from placeholder_vals for each instance of the function.
                kwargs [optional, type=dict]
                    Keyword arguments to pass to the target function.  One or more may have the special value "__placeholder__",
                    which will be replaced with a value from placeholder_vals for each instance of the function.
                max_concurrent [optional, type=int]
                    Maximum number of function instances to run at the same time.  The default is to run an instance for each 
                    value in placeholder_vals at the same time.
    Returns:    A list of the return values from each instance of the function, sorted by the corresponding placeholder value.
    """

    # Dummy function; sets up the pipe for parallelization so the user doesn't have to.
    def doRun(target, pipe, tag, args, kwargs):
        ret_val = target(*args, **kwargs)
        pipe.put((tag, ret_val))
        return

    pipe = Queue(len(list(placeholder_vals)))

    ph_vals = copy.copy(placeholder_vals)

    ret_vals = []

    ph_done = 0

    while len(ph_vals) > 0:
        # We're going to end up popping off chunks of the ph_vals list, so loop until there's nothing left in that list.

        # Do the pop.
        if max_concurrent == -1: max_concurrent = len(ph_vals)
        ph_chunk = [ ph_vals.pop(0) for idx in range(min(max_concurrent, len(ph_vals))) ]

        procs = {}
        for ph_idx, ph_val in enumerate(ph_chunk):
            tag = ph_idx + ph_done
            # Run an instance of the function for every value in this chunk.

            # Replace all the instances of "__placeholder__" in the arguments.
            ph_args = tuple([ ph_val if item == "__placeholder__" else item for item in args ])
            ph_kwargs = dict([ (key, ph_val) if val == "__placeholder__" else (key, val) for key, val in kwargs.iteritems() ])

            # Instantiate the process and start it.  This calls the dummy function doRun above, which calls the actual target function.
            proc = Process(target=doRun, name=str(ph_val), args=(target, pipe, tag, ph_args, ph_kwargs))
            proc.start()
            procs[tag] = proc

        while len(procs) > 0:
            # Wait for the processes to finish, deleting it from the procs list whenever it's finished.  Loops until all the 
            #   processes have finished.
            pipe_out = pipe.get()
            tag, ret_val = pipe_out

            ret_vals.append(pipe_out)
            del procs[tag]

        ph_done += len(ph_chunk)

    # Sort by placeholder value, keep only the return values themselves.
    ret_vals = zip(*sorted(ret_vals, key=lambda x: x[0]))[1]
    if zip_result:
        return zip(*ret_vals)
    else:
        return ret_vals

def probMatchMean(ens, grid_dims=2):
    def PMM(ens, ens_mean):
        ens_mean_ranks = np.argsort(ens_mean.ravel())
        ens_ranks = np.argsort(ens.ravel())
        pm_mean = np.zeros(ens.shape[1:])

        for erm, er in zip(ens_mean_ranks, ens_ranks[::ens.shape[0]]):
            ens_idx = np.unravel_index(er, ens.shape)
            pm_idx = np.unravel_index(erm, pm_mean.shape)

            pm_mean[pm_idx] = ens[ens_idx]
        return pm_mean

    ens_mean = ens.mean(axis=0)

    if len(ens.shape) == 3 or (len(ens.shape) == 4 and grid_dims == 3):
        pm_mean = PMM(ens, ens_mean)
    elif len(ens.shape) == 4 and grid_dims == 2:
        pm_mean = np.empty(ens.shape[1:])
        for wdt in xrange(ens.shape[1]):
            pm_mean[wdt] = PMM(ens[:, wdt], ens_mean[wdt])

    return pm_mean

def publicationFigure(subfigures, layout, corner='ul', colorbar=None):
    import matplotlib
    import pylab
    rows, columns = layout

    base_diag = np.hypot(6, 12)
    size_x, size_y = pylab.gcf().get_size_inches()
    fig_diag = np.hypot(size_x, size_y)

    multiplier = fig_diag / base_diag

    bbox = {'ec':'k', 'fc':'w', 'pad':10, 'lw':multiplier}
    loc = {'ul':(0.0, 1.0), 'ur':(1.0, 1.0)}
    align = {'ul':('left', 'top'), 'ur':('right', 'top')}
    pad_signs = {'l':1, 'r':-1, 'u':-1}

    matplotlib.rcdefaults()
    for rc in ['axes.linewidth', 'xtick.major.size', 'ytick.major.size' ]:
        matplotlib.rcParams[rc] *= multiplier

#   for line in pylab.gca().xaxis.get_ticklines():
#       line.set_linewidth(line.get_linewidth() * multiplier)
#   for line in pylab.gca().yaxis.get_ticklines():
#       line.set_linewidth(line.get_linewidth() * multiplier)

    for idx, sf in enumerate(subfigures):
        pylab.subplot(rows, columns, idx + 1)

        n_row = idx / columns
        n_col = idx % columns
        sf(multiplier=multiplier, layout=(n_row + 1, n_col + 1))

        offset = transforms.ScaledTranslation(pad_signs[corner[1]] * ((bbox['pad'] + bbox['lw']) / (2 * 72.)), pad_signs[corner[0]] * (bbox['pad'] + bbox['lw']) / (2 * 72.), pylab.gcf().dpi_scale_trans)
        text_transform = pylab.gca().transAxes + offset
        text_x, text_y = loc[corner]
        h_align, v_align = align[corner]
        pylab.text(text_x, text_y, "(%s)" % string.ascii_lowercase[idx], transform=text_transform, bbox=bbox, ha=h_align, va=v_align, fontsize=16 * multiplier, fontweight='bold', zorder=1000)

    def onDraw(event):
        min_label_x = np.zeros((columns,), dtype=float)
        for idx in range(len(subfigures)):
            ax = pylab.subplot(rows, columns, idx + 1)
            n_col = idx % columns
            bbox = ax.yaxis.label.get_window_extent().inverse_transformed(ax.transAxes)
            min_label_x[n_col] = min(min_label_x[n_col], bbox.get_points()[:, 0].min())

        for idx in range(len(subfigures)):
            ax = pylab.subplot(rows, columns, idx + 1)
            n_col = idx % columns
            ax.yaxis.label.set_position((min_label_x[n_col], 0.5))

#   pylab.gcf().canvas.mpl_connect('draw_event', onDraw)

    if colorbar:
        bar_label, format, ticks = colorbar[:3]
        tick_labels = colorbar[-1]

        cax = pylab.axes((0.90, 0.1125, 0.020, 0.825))
        bar = pylab.colorbar(cax=cax)
        bar.ax.text(3.5, 0.5, bar_label, rotation=90, transform=bar.ax.transAxes, size=12 * multiplier, va='center')
        bar.set_ticks(ticks)
        bar.set_ticklabels([ format % t for t in tick_labels])
        bar.ax.tick_params(labelsize=12 * multiplier)

    return

if __name__ == "__main__":
    test_array = np.arange(3 * 4 * 5 * 6).reshape((3, 4, 5, 6))
    print test_array.mean(axis=0)
    print probMatchMean(test_array, grid_dims=3)
