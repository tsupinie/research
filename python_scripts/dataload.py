
import Nio as nio

import numpy as np

from util import runConcurrently, decompressVariable
from interp import getInterpFunctions

from computeQuantities import toRecArray, toDict

def _makeZCoordsAGL(z_coords, dz_min=50.):
    z_coords_agl = z_coords - z_coords[1].reshape((1,) + z_coords.shape[1:])
    return z_coords_agl, z_coords

def _buildEnsName(n_ens, t_ens, state='ena'):
    return "%s%03d.hdf%06d" % (state, n_ens + 1, t_ens)

def _buildEnsGrdbas(n_ens, state='ena'):
    return "%s%03d.hdfgrdbas" % (state, n_ens + 1)

def getAxes(base_path, agl=True, z_coord_type=""):
    grdbas_file = _buildEnsGrdbas(0)
    hdf_grdbas = nio.open_file("%s/%s" % (base_path, grdbas_file), mode='r', format='hdf')
    axes = dict([ (ax[:1], decompressVariable(hdf_grdbas.variables[ax])) for ax in ['x', 'y', "zp%s" % z_coord_type] ])
    if agl:
        axes['z'], axes['z_MSL'] = _makeZCoordsAGL(axes['z'])
    else:
        axes['z_AGL'], axes['z'] = _makeZCoordsAGL(axes['z'])

    return axes

def _loadData(hdf_var, bounds):
    if hdf_var[:].max() == np.iinfo(np.int16).max:
        data = decompressVariable(hdf_var, dindex=bounds)
    else:
        data = hdf_var[bounds]
    return data

def loadDataFile(base_path, n_ens, t_ens, derived, interp, wrap=False, prog=True, state='ena', coords='hght', buffer=False):
    var_list, func = derived
    interp_func, bounds_func, axes, points = interp
    if coords == 'hght':
        bounds = slice(*bounds_func(axes, points, coords=coords, buffer=buffer))
        bounds = (bounds, slice(None), slice(None))
        axes['z'] = axes['z'][bounds]

    if coords == 'pres':
        try:
            var_list.remove('p')
        except ValueError:
            pass

        var_list.insert(0, 'p')

    vars = {}

    data_file = _buildEnsName(n_ens, t_ens, state=state)
    hdf = nio.open_file("%s/%s" % (base_path, data_file), mode='r', format='hdf')

    for var in var_list:
        dointerp = False
        if var in hdf.variables and var not in ['x', 'y', 'z']:
            if coords == 'pres' and var == 'p':
                vars[var] = _loadData(hdf.variables[var], slice(None))
                axes['z'] = vars[var]

                bounds = slice(*bounds_func(axes, points, coords=coords))
                bounds = (bounds, slice(None), slice(None))
                vars[var] = vars[var][bounds]
            else:
                vars[var] = _loadData(hdf.variables[var], bounds)

            dointerp = True
        elif var == 'z':
            if 'z_AGL' in axes:
                vars[var] = axes['z_AGL'][bounds]
            elif 'z_MSL' in axes:
                vars[var] = axes['z_MSL'][bounds]

            dointerp = True
        elif var == 'y':
            vars[var] = axes['y'][:]
        elif var == 'x':
            vars[var] = axes['x'][:]
        elif var == 'dx':
            vars[var] = axes['x'][1] - axes['x'][0]
        elif var == 'dy':
            vars[var] = axes['y'][1] - axes['y'][0]

        if points and dointerp and not buffer:
            vars[var] = interp_func(vars[var], axes, points, wrap=wrap)

    result = func(**vars)
    if buffer:
        if result.dtype.fields:
            result = toDict(result)
            for var in result.iterkeys():
                result[var] = interp_func(result[var], axes, points, wrap=wrap)
            result = toRecArray(**result)
        else:
            result = interp_func(result, axes, points, wrap=wrap)
    return result

def _loadEnsembleByTimestep(base_path, members, times, derived, interp, wrap=False, aggregator=None, max_concurrent=-1, state='ena', coords='hght', buffer=False):
    if type(members) in [ int ]:
        members = range(members)
    else:
        members = [ m - 1 for m in members ]

    ensemble = None

    for wdt, t_ens in enumerate(times):
        print "Loading time %d ..." % t_ens
        timestep = runConcurrently(loadDataFile, members, args=(base_path, "__placeholder__", t_ens, derived, interp), kwargs={'wrap':wrap, 'state':state, 'coords':coords, 'buffer':buffer}, max_concurrent=max_concurrent)

        if ensemble is None:
            base = timestep[0]
            ensemble = np.empty((len(members), len(times)) + base.shape, dtype=base.dtype)

        for lde, ens_data in enumerate(timestep):
            ensemble[lde, wdt] = ens_data

    return ensemble

def _loadEnsembleByRun(base_path, members, times, derived, interp, wrap=False, aggregator=None, max_concurrent=-1, state='ena', coords='hght', buffer=False):
    if type(members) in [ int ]:
        members = range(members)
    else:
        members = [ m - 1 for m in members ]

    ensemble = None

    for lde, n_ens in enumerate(members):
        print "Loading member %d ..." % (n_ens + 1)
        run = runConcurrently(loadDataFile, times, args=(base_path, n_ens, "__placeholder__", derived, interp), kwargs={'wrap':wrap, 'state':state, 'coords':coords, 'buffer':buffer}, max_concurrent=max_concurrent)

        if ensemble is None:
            base = run[0]
            ensemble = np.empty((len(members), len(times)) + base.shape, dtype=base.dtype)

        for wdt, ens_data in enumerate(run):
            ensemble[lde, wdt] = ens_data

    return ensemble

def loadEnsemble(base_path, members, times, derived, points=None, wrap=False, agl=False, aggregator=None, max_concurrent=-1, z_coord_type="", coords='hght', fcst=False, buffer=False):
    load_by_run = False
    load_by_run |= type(members) in [ list, tuple ] and len(members) == 1
    load_by_run &= aggregator is None

    states = ['ena', 'enf']

    axes = getAxes(base_path, agl=agl, z_coord_type=z_coord_type)

    interp = getInterpFunctions(axes, points)
    interp = interp + (axes, points)

    if load_by_run: 
        return _loadEnsembleByRun(base_path, members, times, derived, interp, wrap=wrap, max_concurrent=max_concurrent, state=states[int(fcst)])
    else:
        return _loadEnsembleByTimestep(base_path, members, times, derived, interp, wrap=wrap, aggregator=aggregator, max_concurrent=max_concurrent, coords=coords, state=states[int(fcst)], buffer=buffer)

if __name__ == "__main__":
    from computeQuantities import computeReflectivity
#   from util import loadAndInterpolateEnsemble

    from datetime import datetime, timedelta
    import glob
    derived = (['p', 'pt', 'qr', 'qs', 'qh'], computeReflectivity)

    start = datetime.now()
    ens = loadEnsemble("/caps2/tsupinie/1kmf-control/", 40, range(14700, 18300, 300), derived, {'z':1000}, agl=True)
    print "Time to read ensemble (parallel):", datetime.now() - start
    print np.nanmax(ens)

#   files = glob.glob("/caps2/tsupinie/1kmf-control/ena???.hdf014700")
#   files.extend(glob.glob("/caps2/tsupinie/1kmf-control/ena???.hdf01[5678]?00"))
#   start = datetime.now()
#   ens = loadAndInterpolateEnsemble(files, derived[0], derived[1], "/caps2/tsupinie/1kmf-control/ena001.hdfgrdbas", {'z':1000}, agl=True)
#   print "Time to read ensemble (serial):", datetime.now() - start
#   print np.nanmax(ens)
