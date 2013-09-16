
import numpy as np

import scipy.weave as weave

from mpl_toolkits.basemap import Basemap

import Nio as nio

import gc, re, sys
from itertools import izip

# The projection information for all my experiment grids.  These are identical to the arguments that Basemap gets, with one exception:
# ctr_lat and ctr_lon refer to the center of the image (lat_0 and lon_0 refer to the center of the projection, as in the Basemap 
# documentation).

goshen_1km_gs = (1000, 1000)
goshen_1km_proj = {'projection':'lcc', 'resolution':'l',
    'width':(goshen_1km_gs[0] * 255), 'height':(goshen_1km_gs[1] * 255),
    'ctr_lat':41.61975, 'ctr_lon':-104.34843,
    'lat_0':40.61998, 'lon_0':-107.344, 'lat_1':30., 'lon_1':60.
}

goshen_3km_gs = (3000, 3000)
goshen_3km_proj = {'projection':'lcc', 'resolution':'l',
    'width':(goshen_3km_gs[0] * 411), 'height':(goshen_3km_gs[1] * 411),
    'ctr_lat':40.61998, 'ctr_lon':-107.344, 
    'lat_0':40.61998, 'lon_0':-107.344, 'lat_1':30., 'lon_1':60.
}

inflow_stations = {
    14700: {
        'inflow': [ ], 
        'outflow': [ '101A', '102B', '103A', '106B', 'P1', 'P2', 'P3', 'P4', 'P5', 'P7' ], 
        'sounding': [ '88 a', 'Bush', 'NCAR', 'WY/N' ],
    },
    15000: {
        'inflow': [ '104B', '110A', '220B', 'P4' ], 
        'outflow': [ '101A', '102B', '103A', '106B', 'P1', 'P2', 'P3', 'P5' ], 
        'sounding': [ '88 a', 'Bush', 'NCAR', 'WY/N' ],
    },
    15300: {
        'inflow': [ '104B', '105A' ], 
        'outflow': [ '101A', '102B', '103A', '106B', '110A', '214B', '220B', 'P1', 'P2', 'P4' ], 
        'sounding': [ '88 a', 'Bush', 'NCAR', 'WY/N' ],
    },
    15600: {
        'inflow': [ '105A' ], 
        'outflow': [ '101A', '102B', '103A', '104B', '106B', '110A', '214B', '220B', '222B', 'P1', 'P2', 'P4', 'P5' ], 
        'sounding': [ '88 a', 'Bush', 'NCAR', 'WY/N' ],
    },
    15900: {
        'inflow': [ 'P1' ], 
        'outflow': [ '101A', '102B', '103A', '104B', '105A', '106B', '107A', '110A', '214B', '220B', '222B', 'P4', 'P7' ], 
        'sounding': [ '88 a', 'Bush', 'NCAR', 'WY/N' ],
    },
    16200: {
        'inflow': [ '109A', '213A', '216B'], 
        'outflow': [ '101A', '102B', '103A', '104B', '105A', '106B', '107A', '110A', '214B', '220B', '222B', 'P1', 'P4', 'P5', 'P7'  ], 
        'sounding': [ '88 a', 'Bush', 'NCAR', 'WY/N' ],
    },
    16500: {
        'inflow': [ '213A', '216B', '221A' ], 
        'outflow': [ '101A', '102B', '103A', '104B', '105A', '106B', '107A', '109A', '110A', '214B', '217A', '220B', '222B', 'P2', 'P5', 'P7' ], 
        'sounding': [ '88 a', 'Bush', 'NCAR', 'WY/N' ],
    },
    16800: {
        'inflow': [ '213A', '215A', '216B', '217A', '221A' ], 
        'outflow': [ '101A', '102B', '103A', '104B', '105A', '106B', '107A', '109A', '110A', '214B', '220B', '222B', 'P7' ], 
        'sounding': [ '88 a', 'Bush', 'NCAR', 'WY/N' ],
    },
    17100: {
        'inflow': [ '213A', '215A', '221A', 'P2' ], 
        'outflow': [ '101A', '102B', '103A', '104B', '105A', '106B', '107A', '109A', '110A', '214B', '216B', '217A', '220B', '222B', 'P5', 'P7' ], 
        'sounding': [ '88 a', 'Bush', 'NCAR', 'WY/N' ],
    },
    17400: {
        'inflow': [ '215A', 'P2' ], 
        'outflow': [ '101A', '102B', '103A', '104B', '105A', '106B', '107A', '109A', '213A', '214B', '216B', '217A', '220B', '221A', '222B', 'P7' ], 
        'sounding': [ '88 a', 'Bush', 'NCAR', 'WY/N' ],
    },
    17700: {
        'inflow': [ 'P2' ], 
        'outflow': [ '101A', '102B', '103A', '104B', '105A', '106B', '107A', '109A', '213A', '214B', '215A', '216B', '217A', '220B', '221A', '222B', 'P7' ], 
        'sounding': [ '88 a', 'Bush', 'NCAR', 'WY/N' ],
    },
    18000: {
        'inflow': [ 'P2' ], 
        'outflow': [ '101A', '102B', '103A', '104B', '105A', '106B', '107A', '213A', '214B', '215A', '216B', '217A', '220B', '221A', '222B', 'P1', 'P7' ], 
        'sounding': [ '88 a', 'Bush', 'NCAR', 'WY/N' ],
    }
}

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

def _makeZCoordsAGL(z_coords):
    dz_min = 50.
    new_coord_shape = [ 1 ]
    new_coord_shape.extend(z_coords.shape[1:])
    z_coords = z_coords - z_coords[0].reshape(tuple(new_coord_shape))
    return z_coords + dz_min / 2.

def _findZBounds(z_coords, interp_height, decreasing=False):
    argmin_distance = np.argmin(np.abs(z_coords[:, 1:-1, 1:-1] - interp_height), axis=0)
    min_distance = np.empty(argmin_distance.shape)

    for index in np.ndindex(*min_distance.shape):
        new_index = [ argmin_distance[index] ]
        new_index.extend(index)

        min_distance[index] = z_coords[tuple(new_index)]

    if decreasing:
        lower_bound = np.where(min_distance < interp_height, argmin_distance, argmin_distance - 1).min()
        upper_bound = np.where(min_distance > interp_height, argmin_distance, argmin_distance + 1).max() + 1
    else:
        lower_bound = np.where(min_distance > interp_height, argmin_distance, argmin_distance - 1).min()
        upper_bound = np.where(min_distance < interp_height, argmin_distance, argmin_distance + 1).max() + 1

    lower_bound = max(lower_bound, 0)
    upper_bound = min(upper_bound, z_coords.shape[0])

    return lower_bound, upper_bound

def interpolate(data, axes, points, wrap=False):
    global _interp_code
    global _lib_code

    force = 0
    if _interp_code == "":
#       force = 1
        _interp_code = open("interp.c", 'r').read()
#       _interp_code = "#define INTERP_DEBUG\n%s" % _interp_code

    if _lib_code == "":
        _lib_code = open("pylib.c").read()

    data_float = data.astype(np.float32)

    if ('x' in points or 'y' in points) and 'z' in points:
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
    elif 'x' in points or 'y' in points:
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
    elif 'z' in points:
        if type(points['z']) in [int, float]:
            z_coords = axes['z'].astype(np.float32)
            z_pt = points['z']
            data_interp = np.empty(data.shape[1:], dtype=np.float32)
            weave.inline("interpz(data_float_array, z_coords_array, z_pt, data_interp_array, PyObject_IsTrue(wrap) == 1);", 
                ['data_float', 'z_coords', 'z_pt', 'data_interp', 'wrap'], 
                support_code=_lib_code + _interp_code,
                force=force
            )
        else:
            z_coords = axes['z'].astype(np.float32)
            z_pt = points['z'].astype(np.float32)
            data_interp = np.empty(data.shape[1:], dtype=np.float32)
            weave.inline("interpheights(data_float_array, z_coords_array, z_pt_array, data_interp_array, PyObject_IsTrue(wrap) == 1);",
                ['data_float', 'z_coords', 'z_pt', 'data_interp', 'wrap'],
                support_code=_lib_code + _interp_code,
                force=force
            )
    else:
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

def coneHeight(distance, elev_angle, z_base):
    earth_radius = 6371000.
    eff_factor = 4. / 3.
    elev = 3.14592654 / 180. * elev_angle # should be np.pi
    return (np.cos(elev) / np.cos(elev + distance / (earth_radius * eff_factor)) - 1) * earth_radius * eff_factor + z_base;

def computeBounds(z_coords, y_coords, x_coords, points, coords='hght'):
    if points is None:
        lower_bound, upper_bound = 0, z_coords.shape[0]
    else:
        if ('x' in points or 'y' in points) and 'z' in points:
            if coords == 'hght':
                lower_bound, dummy = _findZBounds(z_coords, points['z'].min(), coords=='pres')
                dummy, upper_bound = _findZBounds(z_coords, points['z'].max(), coords=='pres')
            elif coords == 'pres':
                lower_bound, dummy = _findZBounds(z_coords, points['z'].max(), coords=='pres')
                dummy, upper_bound = _findZBounds(z_coords, points['z'].min(), coords=='pres')
        elif 'x' in points or 'y' in points:
            lower_bound, upper_bound = 0, z_coords.shape[0]
        elif 'z' in points:
            lower_bound, upper_bound = _findZBounds(z_coords, points['z'], coords=='pres')
        else:
            max_dist = 0

            for idx in [0, -1]:
                for jdy in [0, -1]:
                    max_dist = max(max_dist, np.hypot(points['x_base'] - x_coords[idx], points['y_base'] - y_coords[jdy]))

            lower_bound, dummy = _findZBounds(z_coords, points['z_base'], coords == 'pres')
            dummy, upper_bound = _findZBounds(z_coords, coneHeight(max_dist, points['elev_angle'], points['z_base']), coords=='pres')
    return lower_bound, upper_bound

def coneHeight(distance, elev_angle, z_base):
    earth_radius = 6371000.
    eff_factor = 4. / 3.
    elev = 3.14592654 / 180. * elev_angle # should be np.pi
    return (np.cos(elev) / np.cos(elev + distance / (earth_radius * eff_factor)) - 1) * earth_radius * eff_factor + z_base;

def computeBounds(z_coords, y_coords, x_coords, points, coords='hght'):
    if points is None:
        lower_bound, upper_bound = 0, z_coords.shape[0]
    else:
        if ('x' in points or 'y' in points) and 'z' in points:
            if coords == 'hght':
                lower_bound, dummy = _findZBounds(z_coords, points['z'].min(), coords=='pres')
                dummy, upper_bound = _findZBounds(z_coords, points['z'].max(), coords=='pres')
            elif coords == 'pres':
                lower_bound, dummy = _findZBounds(z_coords, points['z'].max(), coords=='pres')
                dummy, upper_bound = _findZBounds(z_coords, points['z'].min(), coords=='pres')
        elif 'x' in points or 'y' in points:
            lower_bound, upper_bound = 0, z_coords.shape[0]
        elif 'z' in points:
            lower_bound, upper_bound = _findZBounds(z_coords, points['z'], coords=='pres')
        else:
            max_dist = 0

            for idx in [0, -1]:
                for jdy in [0, -1]:
                    max_dist = max(max_dist, np.hypot(points['x_base'] - x_coords[idx], points['y_base'] - y_coords[jdy]))

            lower_bound, dummy = _findZBounds(z_coords, points['z_base'], coords == 'pres')
            dummy, upper_bound = _findZBounds(z_coords, coneHeight(max_dist, points['elev_angle'], points['z_base']), coords=='pres')
    return lower_bound, upper_bound

def loadGrdbas(grdbas_file, agl):
    grdbas = nio.open_file(grdbas_file, mode='r', format='hdf')
    z_coords = decompressVariable(grdbas.variables['zp'])

    x_coords = decompressVariable(grdbas.variables['x'])
    y_coords = decompressVariable(grdbas.variables['y'])

    if agl:
        z_coords = _makeZCoordsAGL(z_coords)

    return z_coords, y_coords, x_coords

def loadAndInterpolateRun(files, field_list, function, points=None, coord_info=None, grdbas_file=None, agl=True, wrap=False, coords='hght', prog=True):
    times = np.unique([ int(f[-6:]) for f in files ])
    files.sort()

    if coord_info is None:
        if grdbas_file is None:
            print "loadAndInterpolate() needs a grdbas file name or coord information."
            return
        z_coords, y_coords, x_coords, = loadGrdbas(grdbas_file, agl)
        if coords == 'hght':
            lower_bound, upper_bound = computeBounds(z_coords, y_coords, x_coords, points, coords)
    else:
        if coords == 'hght':
            z_coords, y_coords, x_coords, lower_bound, upper_bound = coord_info
        else:
            z_coords, y_coords, x_coords = coord_info

    member_all = 0
    var_ens = {}

    axes = { 'y':y_coords - y_coords[0], 'x':x_coords - x_coords[0] }
    if coords == 'pres':
        if 'p' in field_list:
            field_list.pop(field_list.index('p'))

        field_list.insert(0, 'p')
        field_list.append('z')
    elif coords == 'hght':
        axes['z'] = z_coords[lower_bound:upper_bound]
    else:
        print "Coords '%s' not understood ..." % coords
        sys.exit()

    for wdt, t_ens in enumerate(times):
        if prog:
            print "t_ens = %06d (%s)" % (t_ens, files[wdt])
        hdf = nio.open_file(files[wdt], mode='r', format='hdf')

        if prog:
            sys.stdout.write("Loading ")
            sys.stdout.flush()

        for var in field_list:
            if prog:
                sys.stdout.write("%s ... " % var)
                sys.stdout.flush()

            if var in hdf.variables and var != 'z':
                if coords == 'pres' and var == 'p':
                    pres_coords = hdf.variables[var][:]
                    lower_bound, upper_bound = computeBounds(pres_coords, y_coords, x_coords, points, coords)

                    var_ens[var] = pres_coords[lower_bound:upper_bound]
                    axes['z'] = copy.copy(var_ens[var])
                else:
                    if hdf.variables[var][:].max() == np.iinfo(np.int16).max:
                        bounds = (slice(lower_bound, upper_bound), slice(None), slice(None))
                        var_ens[var] = decompressVariable(hdf.variables[var], dindex=bounds)
                    else:
                        var_ens[var] = hdf.variables[var][lower_bound:upper_bound]

                if points is not None:
                    var_ens[var] = interpolate(var_ens[var], axes, points, wrap=wrap)
            elif var == 'z':
                var_ens[var] = z_coords[lower_bound:upper_bound]
                if points is not None:
                    var_ens[var] = interpolate(var_ens[var], axes, points, wrap=wrap)
            elif var == 'dx':
                var_ens[var] = x_coords[1] - x_coords[0]
            elif var == 'dy':
                var_ens[var] = y_coords[1] - y_coords[0]

            if var in [ 'qr', 'qs', 'qh' ]: var_ens[var] = np.maximum(var_ens[var], np.zeros(var_ens[var].shape))

        if prog:
            sys.stdout.write("\n")
            sys.stdout.flush()

        tmp_member = function(**var_ens)

        if type(member_all) == int:
            member_shape = [ len(times) ]
            member_shape.extend(var_ens[field_list[0]].shape)
            member_all = np.empty(tuple(member_shape), dtype=tmp_member.dtype)

        member_all[wdt] = tmp_member

    return member_all, times

def loadAndInterpolateEnsemble(files, field_list, function, grdbas_file, points=None, agl=True, wrap=False, aggregator=None, coords='hght', prog=True):
    """
    loadAndInterpolateEnsemble()
    Purpose:    Load all members and all times of an ensemble, interpolate them to some space and compute derived quantities.

    Parameters: files [type=list]
                    A list of all the files to load.  The function assumes they have the usual ARPS naming convention and places them in the final array accordingly.
                field_list [type=list]
                    A list of variables to get from all the files (u, v, qv, etc.)
                function [type=function]
                    A function object run on each ensemble member for each time to compute the derived quantity(ies).  The only arguments to the function should be 
                        **kwargs.  The data type of the final array will take the form of the data type of the return value from this function, so you can return
                        multiple derived quantities using a numpy record array.
                grdbas_file [type=str]
                    The name of the grid-base file to use when looking for the grid coordinates.
                points [type=dict, optional]
                    The points to interpolate to.  This should be a dictionary and can take several forms.  For point interpolations, it should contain x, y, and z
                        keys, the values for which should be the x, y, and z coordinates for each point, respectively.  The 'z' key should point to pressures in the 
                        case that coords is given as 'pres'.  If this argument is not specified, no interpolation will be done.  
                agl [type=boolean, optional]
                    Whether the z-coordinates in the points should be interpreted as AGL (True) or MSL (False).  AGL (True) is the default.
                wrap [type=boolean, optional]
                    Whether or not to "wrap" out-of-bounds values in the interpolation.  If True, will assign the value of the closest point in the domain.  If False,
                        assigns NaN.
                aggregator [type=function, optional]
                    If given, aggregates all ensemble members using the specified funtion (e.g. np.mean()) after interpolation and before computing the derived quantity.
                coords [type=str, optional]
                    The height coordinates to use for the interpolation.  Default is 'hght', which will use height coordinates from the grid-base file.  Could also be
                        'pres', which means to use pressure from the history file for the vertical coordinates.  In this case, the 'z' key in the points dictionary should
                        point to pressures of the points and not heights.
                prog [type=boolean, optional]
                    Specifies whether to show the progress for loading variables from the file.  Default is True.  Intended for the case where you're running multiple 
                        instances of this function at the same time, and you don't want the progress readouts stepping on each other.

    Returns:    In the case where aggregator is not specified and points is not specified, it returns a numpy array with dimensions NE x NT x NZ x NY x NX, where NE is
                    the number of ensemble members, NT is the number of time steps, NZ is the number of z points, NY is the number of y points, and NX is the number of x
                    points.  
                If aggregator is specified, the function will return two arrays.  The first array will not have an NE dimension, and will be the aggregated ensemble.  The
                    second array will be the result of computing the derived quantity (see the function argument) on all ensemble members, and will be NE x NT x NZ x NY x
                    NX  
                If points is specified, the resulting array will be depend on the interpolation method used.  For point interpolation, it will be NE x NT x NP, where 
                    NP is the number of points.
                If the two are given in combination, it *should* be smart enough to handle this correctly.  In that case, the dimensions of the output array will be
                    some combination of the above cases.
                For all cases of input parameters, the list of ensemble members and time steps are returned, as well.
    """
    files_partitioned = []
    ens_members = np.unique([ f[-13:-10] for f in files ])

    for n_ens in ens_members:
        files_partitioned.append([ f for f in files if re.search(r"e[\w]{2}%s" % n_ens, f) ])

    n_ensemble_members = len(ens_members)

    z_coords, y_coords, x_coords = loadGrdbas(grdbas_file, agl)
    if coords == 'hght':
        lower_bound, upper_bound = computeBounds(z_coords, y_coords, x_coords, points, coords)
        coord_info = (z_coords, y_coords, x_coords, lower_bound, upper_bound)
    else:
        coord_info = (z_coords, y_coords, x_coords)

    param_all = 0
    ens_times = None

    if aggregator is None:
        quantity_func = function
    else:
        quantity_func = lambda **d: np.array(zip(*[ d[k].ravel() for k in sorted(d.keys()) ]), dtype=[ (k, float) for k in sorted(d.keys()) ]).reshape(d[d.keys()[0]].shape)

    for lde, n_ens in enumerate(ens_members):
        if prog:
            print "n_ens = %s" % n_ens
        param, ens_times = loadAndInterpolateRun(files_partitioned[lde], copy.copy(field_list), quantity_func, points, coord_info=coord_info, wrap=wrap, coords=coords, prog=prog)

        if type(param_all) == int:
            param_shape = [ len(ens_members) ]
            param_shape.extend(param.shape)
            param_all = np.empty(tuple(param_shape), dtype=param.dtype)

        param_all[lde] = param
        gc.collect()

    if aggregator is not None:
        aggregate_dict = {}
        param_dict = {}
        for var in param_all.dtype.fields.iterkeys():
            aggregate_dict[var] = aggregator(param_all[var])
            param_dict[var] = param_all[var]

        agg_param_all = function(**aggregate_dict)
        param_all = function(**param_dict)
        return agg_param_all, param_all, ens_members, ens_times
    else:
        return param_all, ens_members, ens_times

def setupMapProjection(proj, grid_spacing, bounds=None):
    if bounds is None:
        bound_x, bound_y = slice(None), slice(None)
    else:
        bound_x, bound_y = bounds

    new_proj = copy.copy(proj)

    dx, dy = grid_spacing
    nx = new_proj['width'] / dx
    ny = new_proj['height'] / dy

    if 'ctr_lat' in new_proj:
        ctr_lat = new_proj['ctr_lat']
        del new_proj['ctr_lat']

    if 'ctr_lon' in new_proj:
        ctr_lon = new_proj['ctr_lon']
        del new_proj['ctr_lon']

    map = Basemap(**new_proj) 

    ctr_x, ctr_y = map(ctr_lon, ctr_lat)
    move_x = ctr_x - new_proj['width'] / 2
    move_y = ctr_y - new_proj['height'] / 2

    lower_bound_x = bound_x.start
    if lower_bound_x is None: lower_bound_x = 0
    lower_bound_y = bound_y.start
    if lower_bound_y is None: lower_bound_y = 0
    upper_bound_x = bound_x.stop
    if upper_bound_x is None: upper_bound_x = nx
    upper_bound_y = bound_y.stop
    if upper_bound_y is None: upper_bound_y = ny

    llcrnrlon, llcrnrlat = map(move_x + dx * lower_bound_x, move_y + dy * lower_bound_y, inverse=True)
    urcrnrlon, urcrnrlat = map(move_x + dx * upper_bound_x, move_y + dy * upper_bound_y, inverse=True)

    del new_proj['width']
    del new_proj['height']

    new_proj['llcrnrlat'] = llcrnrlat
    new_proj['llcrnrlon'] = llcrnrlon
    new_proj['urcrnrlat'] = urcrnrlat
    new_proj['urcrnrlon'] = urcrnrlon

    return new_proj

def plot_map(plot_data, proj_info, title, file_name=None, color_bar='refl', topo=None, vectors=None, obs=None, pcolormesh=False, subplot=None):
    import matplotlib
    matplotlib.use('agg')
    import pylab

    if subplot is None:
        pylab.clf()
    else:
        pylab.subplot(subplot)
    nx, ny = plot_data.shape

    map = Basemap(**proj_info) 

    llcrnrx, llcrnry = map(proj_info['llcrnrlon'], proj_info['llcrnrlat'])
    urcrnrx, urcrnry = map(proj_info['urcrnrlon'], proj_info['urcrnrlat'])

    dx = (urcrnrx - llcrnrx) / nx
    dy = (urcrnry - llcrnry) / ny

    if topo is not None:
        topo_data, topo_lats, topo_lons = topo

        topo_lats, topo_lons = np.meshgrid(topo_lats, topo_lons)
        topo_x, topo_y = map(topo_lons, topo_lats)

        map.contourf(topo_x, topo_y, topo_data, cmap=pylab.get_cmap('gray'))

    if color_bar == 'refl':
        levels = range(10, 80, 10)
        color_map = pylab.get_cmap('jet')
    elif color_bar == 'drefl':
        levels = range(-30, 35, 5)
        color_map = pylab.get_cmap('RdBu')
    elif color_bar == 'radv':
        levels = range(-30, 30, 5)
        color_map = pylab.get_cmap('RdBu')
    elif color_bar == 'pt':
        levels = range(296, 321, 2)
        color_map = pylab.get_cmap('jet')
    elif color_bar == 'prob':
        levels = np.arange(0.1, 1.1, 0.1)
        color_map = pylab.get_cmap('RdYlBu_r')
    elif color_bar == 'time':
        levels = np.unique(plot_data.ravel())[1:]
        color_map = pylab.get_cmap('Accent')

    x, y = np.meshgrid(dx * np.arange(nx), dy * np.arange(ny))
    if pcolormesh:
        color_map.set_under('#ffffff')
        pylab.pcolormesh(x, y, plot_data, cmap=color_map, vmin=levels[0], vmax=levels[-1])
    else:
        map.contourf(x, y, plot_data, levels=levels, cmap=color_map)

    pylab.colorbar()

    map.drawcoastlines(linewidth=1.0)
    map.drawstates(linewidth=1.0)
    map.readshapefile("countyp020", "counties", linewidth=0.5)

    if obs is not None:
        for stid, ob in obs.iteritems():
            potential_temperature = (5. / 9. * (ob['temperature'] - 32) + 273.15) * (29.5250192 / ob['pressure']) ** (2. / 7.)
            ob_x, ob_y = map(ob['Longitude'], ob['Latitude'])

            ob_ax_x, ob_ax_y = (pylab.gca().transData + pylab.gca().transAxes.inverted()).transform(np.array([ob_x, ob_y]))

            if ob_ax_x > 0 and ob_ax_x <= 1 and ob_ax_y > 0 and ob_ax_y <= 1:
                pylab.gca().add_patch(Circle((ob_x, ob_y), 4000, fc='none', ec='k'))
                pylab.text(ob_x, ob_y, "%5.1f" % potential_temperature, size='x-small', ha='right', va='bottom')

    if vectors is not None:
        stride = 15
        u, v = vectors
        pylab.quiver(x[::stride, ::stride], y[::stride, ::stride], u[::stride, ::stride], v[::stride, ::stride])

    pylab.title(title)

    if file_name is not None:
        pylab.savefig(file_name)
    return

def drawPolitical(map, scale_len=0):
    """
    drawPolitical()
    Purpose:    Draws all the political boundaries on a map.  Has the additional advantage that it doesn't require loading in all the shapefiles
                    if you want to draw the boundaries again.  Also makes use of a states shapefile that's not internal to Basemap because the 
                    one in basemap is off by a noticeable amount.  If you want it, ask me.
    Parameters: map [type=Basemap]
                    A basemap object to use when drawing boundaries and loading shapefiles.
                scale_len [type=int, optional]
                    Specifies a length in km for a scale marker drawn in the lower-right corner of the map.  Default is to not draw a marker.
    Returns:    [nothing]
    """
    import pylab
    from matplotlib.patches import Polygon
    map.drawcountries(linewidth=1.0)
    map.drawcoastlines(linewidth=1.0)

    state_list = ['CO', 'NE', 'WY', 'SD', 'MT', 'ID', 'UT', 'NV', 'AZ', 'NM', 'TX', 'OK', 'KS' ]

    if not hasattr(map, 'counties'):
        map.readshapefile("countyp020", 'counties', linewidth=0.5)
    else:
        for county, data in izip(map.counties_info, map.counties):
            if county['STATE'] in state_list:
                pylab.gca().add_patch(Polygon(data, ec='k', fc='none', linewidth=0.5))

    if not hasattr(map, 'states'):
        map.readshapefile("tl_2011_us_state", 'states', linewidth=1.5)
    else:
        for state, data in izip(map.states_info, map.states):
            if state['STUSPS'] in state_list:
                pylab.gca().add_patch(Polygon(data, ec='k', fc='none', linewidth=1.5))

    if scale_len > 0:
        pylab.plot([scale_len * 1000, scale_len * 2000], [scale_len * 1000, scale_len * 1000], 'k-', lw=2.5)
        pylab.text(scale_len * 1500, scale_len * 800, "%d km" % scale_len, ha='center', va='top')
    return

print "*** WARNING ***"
print "This script is using legacy code.  Please re-write to use the grid and temporal objects, dataload functions, etc."
