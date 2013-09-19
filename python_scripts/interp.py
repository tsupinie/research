
import numpy as np
import scipy.weave as weave

_weave_path = "/data6/tsupinie/weave/"

_interp_code = ""
_lib_code = ""

def _findZBounds(z_coords, interp_height, decreasing=False, buffer=False):
    kdz_buffer = 0
    if buffer: kdz_buffer = 1

    argmin_distance = np.argmin(np.abs(z_coords[:, 1:-1, 1:-1] - interp_height), axis=0)
    min_distance = np.empty(argmin_distance.shape)

    for index in np.ndindex(*min_distance.shape):
        min_distance[index] = z_coords[(argmin_distance[index],) + index]

    if decreasing:
        lower_bound = np.where(min_distance > interp_height, argmin_distance, argmin_distance - 1).min()
        upper_bound = np.where(min_distance < interp_height, argmin_distance, argmin_distance + 1).max() + 1
    else:
        lower_bound = np.where(min_distance < interp_height, argmin_distance, argmin_distance - 1).min()
        upper_bound = np.where(min_distance > interp_height, argmin_distance, argmin_distance + 1).max() + 1

    lower_bound = max(lower_bound, 0)
    upper_bound = min(upper_bound, z_coords.shape[0])

    return lower_bound - 1, upper_bound + 1

def coneHeight(distance, elev_angle, z_base):
    earth_radius = 6371000.
    eff_factor = 4. / 3.
    elev = 3.14592654 / 180. * elev_angle # should be np.pi
    return (np.cos(elev) / np.cos(elev + distance / (earth_radius * eff_factor)) - 1) * earth_radius * eff_factor + z_base;

#
# Code for interpolating to points ...
#
def _interpolatePoints(data, axes, points, wrap=False, force=False):
    data_float = data.astype(np.float32)

    z_coords = axes['z'].astype(np.float32)
    y_coords = axes['y'].astype(np.float32)
    x_coords = axes['x'].astype(np.float32)

    z_pts = points['z'].astype(np.float32)
    y_pts = points['y'].astype(np.float32)
    x_pts = points['x'].astype(np.float32)

    data_interp = np.empty(z_pts.shape, dtype=np.float32)
    weave.inline("interppts(data_float_array, z_coords_array, z_pts_array, y_coords_array, y_pts_array, x_coords_array, x_pts_array, data_interp_array, PyObject_IsTrue(wrap) == 1);", 
        ['data_float', 'z_coords', 'z_pts', 'y_coords', 'y_pts', 'x_coords', 'x_pts', 'data_interp', 'wrap'], 
        support_code=_lib_code + _interp_code,
        force=force
    )

    return data_interp

def _computeBoundsPoints(axes, points, coords='hght', buffer=False):
    lower_bound, dummy = _findZBounds(axes['z'], points['z'].min(), coords=='pres', buffer=buffer)
    dummy, upper_bound = _findZBounds(axes['z'], points['z'].max(), coords=='pres', buffer=buffer)
    return lower_bound, upper_bound
#
# Code for interpolating to columns ...
#
def _interpolateColumn(data, axes, points, wrap=False, force=False):
    data_float = data.astype(np.float32)

    y_coords = axes['y'].astype(np.float32)
    x_coords = axes['x'].astype(np.float32)

    y_pt = np.float32(points['y'])
    x_pt = np.float32(points['x'])

    data_interp = np.empty(data.shape[0], dtype=np.float32)

    weave.inline("interpsounding(data_float_array, y_coords_array, y_pt, x_coords_array, x_pt, data_interp_array);",
        ['data_float', 'y_coords', 'y_pt', 'x_coords', 'x_pt', 'data_interp'],
        support_code=_lib_code + _interp_code,
        force=force
    )

    return data_interp

def _computeBoundsColumn(axes, points, coords='hght', buffer=False):
    return 0, axes['z'].shape[0]

#
# Code for interpolating to a single level ...
#
def _interpolateLevel(data, axes, points, wrap=False, force=False):
    data_float = data.astype(np.float32)

    z_coords = axes['z'].astype(np.float32)
    z_pt = points['z']
    data_interp = np.empty(data.shape[1:], dtype=np.float32)
    weave.inline("interpz(data_float_array, z_coords_array, z_pt, data_interp_array, PyObject_IsTrue(wrap) == 1);", 
        ['data_float', 'z_coords', 'z_pt', 'data_interp', 'wrap'], 
        support_code=_lib_code + _interp_code,
        force=force
    )
    return data_interp

def _computeBoundsLevel(axes, points, coords='hght', buffer=False):
    return _findZBounds(axes['z'], points['z'], coords=='pres', buffer=buffer)

#
# Code for interpolating to a different height for each column ...
#
def _interpolateHeights(data, axes, points, wrap=False, force=False):
    data_float = data.astype(np.float32)

    z_coords = axes['z'].astype(np.float32)
    z_pt = points['z'].astype(np.float32)
    data_interp = np.empty(data.shape[1:], dtype=np.float32)
    weave.inline("interpheights(data_float_array, z_coords_array, z_pt_array, data_interp_array, PyObject_IsTrue(wrap) == 1);",
        ['data_float', 'z_coords', 'z_pt', 'data_interp', 'wrap'],
        support_code=_lib_code + _interp_code,
        force=force
    )
    return data_interp

def _computeBoundsHeights(axes, points, coords='hght', buffer=False):
    lower_bound, dummy = _findZBounds(axes['z'], points['z'].min(), coords=='pres', buffer=buffer)
    dummy, upper_bound = _findZBounds(axes['z'], points['z'].max(), coords=='pres', buffer=buffer)
    return lower_bound, upper_bound

#
# Code for interpolating to a cone ...
#
def _interpolateCone(data, axes, points, wrap=False, force=False):
    data_float = data.astype(np.float32)

    z_coords = axes['z'].astype(np.float32)
    y_coords = axes['y'].astype(np.float32)
    x_coords = axes['x'].astype(np.float32)
       
    z_base = np.float32(points['z_base'])
    y_base = np.float32(points['y_base'])
    x_base = np.float32(points['x_base'])
    elev_angle = np.float32(points['elev_angle'])

    data_interp = np.empty(data.shape[1:], dtype=np.float32)
    beam_height = np.empty(data.shape[1:], dtype=np.float32)

    weave.inline("interpcone(data_float_array, z_coords_array, y_coords_array, x_coords_array, elev_angle, z_base, y_base, x_base, data_interp_array, beam_height_array, PyObject_IsTrue(wrap) == 1);",
        ['data_float', 'z_coords', 'y_coords', 'x_coords', 'z_base', 'y_base', 'x_base', 'elev_angle', 'data_interp', 'beam_height', 'wrap'],
        support_code=_lib_code + _interp_code,
        force=force,
        libraries=['m'],
    )
    return data_interp

def _computeBoundsCone(axes, points, coords='hght', buffer=False):
    max_dist = 0

    for idx in [0, -1]:
        for jdy in [0, -1]:
            max_dist = max(max_dist, np.hypot(points['x_base'] - axes['x'][idx], points['y_base'] - axes['y'][jdy]))

    lower_bound, dummy = _findZBounds(axes['z'], points['z_base'], coords == 'pres', buffer=buffer),
    dummy, upper_bound = _findZBounds(axes['z'], coneHeight(max_dist, points['elev_angle'], points['z_base']), coords=='pres', buffer=buffer)

    return lower_bound, upper_bound  

def _interpolateLayer(data, axes, points, wrap=False, force=False):
    data_float = data.astype(np.float32)
    lbound, ubound = _computeBoundsLayer(axes, points)

    z_coords = axes['z'].astype(np.float32)
    z_lb, z_ub = points['z']
    data_interp = np.empty((ubound - lbound,) + data.shape[1:], dtype=np.float32)
    z_interp = np.empty((ubound - lbound,) + data.shape[1:], dtype=np.float32)
    weave.inline("interplayer(data_float_array, z_coords_array, z_lb, z_ub, data_interp_array, z_interp_array, PyObject_IsTrue(wrap) == 1);", 
        ['data_float', 'z_coords', 'z_lb', 'z_ub', 'data_interp', 'z_interp', 'wrap'], 
        support_code=_lib_code + _interp_code,
        force=1
    )
    return data_interp

def _computeBoundsLayer(axes, points, coords='hght', buffer=False):
    lb, ub = points['z']
    lower_bound, dummy = _findZBounds(axes['z'], lb, coords=='pres', buffer=buffer)
    dummy, upper_bound = _findZBounds(axes['z'], ub, coords=='pres', buffer=buffer)
    return lower_bound, upper_bound

#
# Code for no interpolations
#
def _interpolateNone(data, axes, points, wrap=False, force=False):
    return data

def _computeBoundsNone(axes, points, coords='hght', buffer=False):
    return 0, axes['z'].shape[0]

def _interpolateSigma(data, axes, points, wrap=False, force=False):
    return data[0]

def _computeBoundsSigma(axes, points, coords='hght', buffer=False):
    kdz_buffer = 0
    if buffer: kdz_buffer = 1
    model_level = points['sigma']
    return model_level - kdz_buffer, model_level + kdz_buffer + 1

def getInterpFunctions(axes, points, wrap=False, debug=False):
    global _interp_code
    global _lib_code

    if _interp_code == "":
        _interp_code = open("%s/interp.c" % _weave_path, 'r').read()
        if debug:
            _interp_code = "#define INTERP_DEBUG\n%s" % _interp_code

    if _lib_code == "":
        _lib_code = open("%s/pylib.c" % _weave_path, 'r').read()

    if points is None:
        return _interpolateNone, _computeBoundsNone
    elif ('x' in points or 'y' in points) and 'z' in points:
        return _interpolatePoints, _computeBoundsPoints 
    elif 'x' in points or 'y' in points:
        return _interpolateColumn, _computeBoundsColumn
    elif 'z' in points:
        if type(points['z']) in [int, float]:
            return _interpolateLevel, _computeBoundsLevel
        elif type(points['z']) in [ tuple ]:
            return _interpolateLayer, _computeBoundsLayer
        else:
            return _interpolateHeights, _computeBoundsHeights
    elif 'sigma' in points:
        return _interpolateSigma, _computeBoundsSigma
    else:
        return _interpolateCone, _getBoundsCone
