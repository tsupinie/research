
import numpy as np

from arpsmodelobs import ARPSModelObsFile
from radarobsfile import RadarObsFile
from grid import goshen_1km_grid
from util import publicationFigure

from collections import OrderedDict

import matplotlib
#matplotlib.use('agg')
import pylab

def main():
    base_path = "/caps2/tsupinie/"
    experiments = OrderedDict([('1kmf-sndr0h=25km', 'CTRL'), ('1kmf-zs25-offtime-05XP', 'CTRL_OFF')])
    times = {'1kmf-sndr0h=25km':'220000', '1kmf-zs25-offtime-05XP':'215719'}

    bounds = (slice(110, 135), slice(118, 143))
    grid = goshen_1km_grid(bounds)    
    bounds = grid.getBounds()

    radar_obs = []
    model_obs = []
    for exp in experiments.iterkeys():
        radar_obs.append(RadarObsFile("qc/1km/vols/05XP.20090605.%s" % times[exp]))
        try:
            mo = ARPSModelObsFile("%s/%s/05XPbg014400" % (base_path, exp))
        except AssertionError:
            mo = ARPSModelObsFile("%s/%s/05XPbg014400" % (base_path, exp), mpi_config=(2, 12))
        except:
            print "Couldn't load radial velocity."
            mo = {'vr':np.zeros((1, 255, 255), dtype=np.float32)}

        model_obs.append(mo)

    xs, ys = grid.getXY()
    def subplotFactory(radar_obs, model_obs):
        vr_obs = np.where(radar_obs['Z'] < 15, np.nan, radar_obs['vr'])
        def doSubplot(multiplier=1.0, layout=(-1, -1)):
            pylab.contour(xs, ys, vr_obs[0][bounds], colors='k', linestyles='-', levels=np.arange(-50, 60, 10))
            pylab.contour(xs, ys, model_obs['vr'][0][bounds], colors='k', linestyles='--', levels=np.arange(-50, 60, 10))
            pylab.contourf(xs, ys, (model_obs['vr'][0] - vr_obs[0])[bounds], cmap=matplotlib.cm.get_cmap('RdBu_r'), zorder=-1)

            grid.drawPolitical()
            return
        return doSubplot 

    subplots = []
    for robs, mobs in zip(radar_obs, model_obs):
        subplots.append(subplotFactory(robs, mobs))

    pylab.figure(figsize=(12, 8))
    publicationFigure(subplots, (1, 2), corner='ur', colorbar=("$v_r$ difference (m s$^{-1}$)", "%d", np.arange(-10, 11, 1)))
    pylab.savefig("radar_obs_comparison.png")
    pylab.close()
    return

if __name__ == "__main__":
    main()
