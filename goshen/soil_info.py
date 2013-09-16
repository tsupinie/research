
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

from datetime import datetime, timedelta

from grid import goshen_1km_grid, goshen_3km_grid
from dataload import loadEnsemble
from computeQuantities import toRecArray
from arpsmodelobs import ARPSModelObsFile

def plotSoil(ens_mean, ens_spread, grid, title, file_name, refl=None):
    xs, ys = grid.getXY()

    print ens_mean.shape

    pylab.contourf(xs, ys, ens_mean, cmap=matplotlib.cm.get_cmap('jet'))#, levels=np.arange(0.001, 0.015, 0.001))
    pylab.colorbar()
    CS = pylab.contour(xs, ys, ens_spread, colors='k')#, levels=np.arange(0.09, 0.36, 0.03))
    pylab.clabel(CS, fontsize='x-small', inline_spacing=0)
    if refl is not None:
        pylab.contour(xs, ys, refl, colors='k', levels=np.arange(20, 80, 20))

    grid.drawPolitical()
    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def main():
    base_path = "/caps1/tsupinie/"
    n_ens_members = [ 21 ]
    t_ens_start = 0
    t_ens_stop = 0
    t_ens_step = 300

    times = range(t_ens_start, t_ens_stop + t_ens_step, t_ens_step)
    ens = loadEnsemble(base_path, n_ens_members, times, (['tsoil', 'qsoil'], toRecArray), z_coord_type="soil")
    ens_qv = loadEnsemble(base_path, n_ens_members, times, (['qv'], lambda **e: e['qv']), points={'z':75}, agl=True)

    ens_spread = np.empty(ens.shape[1:], dtype=ens.dtype)
    ens_mean = np.empty(ens.shape[1:], dtype=ens.dtype)
    ens_qv_mean = ens_qv.mean(axis=0)

    for field in ens.dtype.fields.keys():
        ens_spread[field] = ens[field].std(axis=0, ddof=1)
        ens_mean[field] = ens[field].mean(axis=0)
    
    base_time = datetime(2009, 6, 5, 18, 0, 0)
    grid = goshen_3km_grid()

    for wdt, time in enumerate(times):
#       try:
#           mo = ARPSModelObsFile("%s/KCYSan%06d" % (base_path, time))
#       except AssertionError:
#           mo = ARPSModelObsFile("%s/KCYSan%06d" % (base_path, time), mpi_config=(2, 12))

        dt = base_time + timedelta(seconds=time)
        plotSoil(ens_mean['tsoil'][wdt, 0, 1], ens_qv_mean[wdt], grid, "Deep Soil Temperature at %s UTC (nstyp=0)" % dt.strftime("%H%M"), "tsoil_%06d_nstyp=0_deep.png" % time)#, refl=mo['Z'][0])
        plotSoil(ens_mean['tsoil'][wdt, 1, 1], ens_qv_mean[wdt], grid, "Deep Soil Temperature at %s UTC (nstyp=1)" % dt.strftime("%H%M"), "tsoil_%06d_nstyp=1_deep.png" % time)#, refl=mo['Z'][0])
        plotSoil(ens_mean['tsoil'][wdt, 0, 0], ens_qv_mean[wdt], grid, "Surface Soil Temperature at %s UTC (nstyp=0)" % dt.strftime("%H%M"), "tsoil_%06d_nstyp=0_sfc.png" % time)#, refl=mo['Z'][0])
        plotSoil(ens_mean['tsoil'][wdt, 1, 0], ens_qv_mean[wdt], grid, "Surface Soil Temperature at %s UTC (nstyp=1)" % dt.strftime("%H%M"), "tsoil_%06d_nstyp=1_sfc.png" % time)#, refl=mo['Z'][0])
    return

if __name__ == "__main__":
    main()
