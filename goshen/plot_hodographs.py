
import numpy as np
from scipy.interpolate import interp1d

from dataload import getAxes
from temporal import goshen_1km_temporal
from grid import goshen_1km_grid
from util import publicationFigure
from arpsmodelobs import ARPSModelObsFile
from color_tables import NWSRef

import matplotlib
matplotlib.use('agg')
import pylab
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import cPickle
from collections import OrderedDict

def plotHodoBackground():
    pylab.gca().set_aspect('equal', adjustable='box')
    pylab.gca().axes.get_xaxis().set_ticks([])
    pylab.gca().axes.get_yaxis().set_ticks([])

    for w_spd in xrange(5, 60, 5):
        ring = Circle((0, 0), w_spd, edgecolor='#444444', facecolor='none', linestyle='dashed', zorder=0)
        pylab.gca().add_patch(ring)

        if not (w_spd % 10):
            pylab.text(w_spd, -1, "%d" % w_spd, color='#444444', ha='center', va='top', clip_on=True, clip_box=pylab.gca().get_position(), zorder=0, size=22)

    pylab.axvline(x=0, color='k', zorder=0)
    pylab.axhline(y=0, color='k', zorder=0)

    pylab.xlim((-15, 35))
    pylab.ylim((-15, 35))
    return

def main():
    base_path = "/caps2/tsupinie/"
    min_ens = 17
    experiments = OrderedDict([('1kmf-sndr0h=25km', 'CTRL'), ('1kmf-zs25-no-05XP', 'NO_MWR'), ('1kmf-z-no-snd', 'NO_SND'), ('1kmf-zs25-no-mm', 'NO_MM'), ('1kmf-zs25-no-mm-05XP', 'NO_MWR_MM'), ('1kmf-z-no-v2', 'NO_V2')])

    temp = goshen_1km_temporal(start=14400, end=18000)
    bounds = (slice(120, 195), slice(95, 170))
    grid = goshen_1km_grid(bounds)
    bounds = grid.getBounds()

    obs = cPickle.load(open('soundings.pkl', 'r'))[2]
    u_obs, v_obs, alt = obs['u_wind'], obs['v_wind'], obs['altitude']
    good_uv = np.where((u_obs != 9999.0) & (v_obs != 9999.0) & (alt != 99999.0))
    u_obs, v_obs, alt = u_obs[good_uv], v_obs[good_uv], alt[good_uv]

    hodo_data = []
    z_axes = []
    for exp in experiments.iterkeys():
        hodo_data.append(cPickle.load(open("hodo_pkl/%s_hodo.pkl" % exp, 'r')))
        z_axes.append(getAxes("%s%s" % (base_path, exp), agl=False)['z'][:, 105, 180])

    def subplotFactory(exp, hodo_data, z_axis):
        def doSubplot(multiplier=1.0, layout=(-1, -1)):
            plt_points = np.where(z_axis < 10000.)[0]

            for n_ens in xrange(hodo_data.shape[0]):
                pylab.plot(hodo_data['u'][n_ens, plt_points], hodo_data['v'][n_ens, plt_points], color='b', zorder=-1)

            pylab.plot(hodo_data['u'][:, plt_points].mean(axis=0), hodo_data['v'][:, plt_points].mean(axis=0), color='k', lw=1.5, zorder=-1)
            plotHodoBackground()
            return
        return doSubplot

#   for wdt, time_sec in enumerate(temp):
#       subplots = []
#       for exp, hodo, z_axis in zip(experiments.iterkeys(), hodo_data, z_axes):
#           subplots.append(subplotFactory(exp, hodo[:, wdt], z_axis))
#
#       pylab.figure(figsize=(12, 8))
#       pylab.subplots_adjust(left=0.05, bottom=0.1, right=0.875, top=0.975, hspace=0.1, wspace=0.1)
#       publicationFigure(subplots, (2, 3), corner='ur')
#       pylab.savefig("hodographs_%d.png" % time_sec)

    xs, ys = grid.getXY()
    colors = dict(zip(experiments.keys(), ['k', 'r', 'g', 'b', 'c', 'm']))
    for wdt, time_sec in enumerate(temp):
        try:
            mo = ARPSModelObsFile("%s/%s/KCYS%03dan%06d" % (base_path, experiments.keys()[0], min_ens, time_sec))
        except AssertionError:
            mo = ARPSModelObsFile("%s/%s/KCYS%03dan%06d" % (base_path, experiments.keys()[0], min_ens, time_sec), mpi_config=(2, 12))
        except:
            print "Can't load reflectivity ..."
            mo = {'Z':np.zeros((1, 255, 255), dtype=np.float32)}

        pylab.figure(figsize=(12, 12))

        plt_obs = np.where((alt - alt[0] >= 50) & (alt < 10000.))[0]
        obs_u_marker = interp1d(alt, u_obs)(1000 * np.arange(2, 10, 2))
        obs_v_marker = interp1d(alt, v_obs)(1000 * np.arange(2, 10, 2))

        pylab.plot(u_obs[plt_obs], v_obs[plt_obs], color='#333333', lw=3, ls='--', label='Observed')
        pylab.plot(u_obs[plt_obs], v_obs[plt_obs], color='#333333', lw=1, ls='-')
        pylab.plot(obs_u_marker, obs_v_marker, marker='o', color="#333333", linestyle='none', markersize=12)

        for exp, hodo, z_axis in zip(experiments.iterkeys(), hodo_data, z_axes):
            u_mean = hodo['u'][:, wdt].mean(axis=0)
            v_mean = hodo['v'][:, wdt].mean(axis=0)

            u_marker = interp1d(z_axis, u_mean)(1000 * np.arange(2, 10, 2))
            v_marker = interp1d(z_axis, v_mean)(1000 * np.arange(2, 10, 2))

            plt_points = np.where((z_axis - z_axis[0] >= 50) & (z_axis < 10000.))[0]
            pylab.plot(u_mean[plt_points], v_mean[plt_points], color=colors[exp], label=experiments[exp], lw=2, zorder=-1)
            pylab.plot(u_marker, v_marker, "%so" % colors[exp], ms=12)

        plotHodoBackground()

        pylab.legend(loc=2, prop={'size':18})
        pylab.xlabel('u (m s$^{-1}$)', size=18)
        pylab.ylabel('v (m s$^{-1}$)', size=18)

        par_ax = pylab.gca()
        ax_width_inches, ax_height_inches = [ f * a for f, a in zip(pylab.gcf().get_size_inches(), par_ax.get_position().size) ]
        ins_ax = inset_axes(pylab.gca(), 0.3 * ax_width_inches, 0.3 * ax_width_inches, loc=1)
        pylab.contourf(xs, ys, mo['Z'][0][bounds], cmap=NWSRef, levels=np.arange(5, 80, 5))

        pylab.plot(1000 * (164 - bounds[1].start), 1000 * (103 - bounds[0].start), 'k*', ms=24)

        grid.drawPolitical()

        pylab.savefig("hodographs_comb_%d.png" % time_sec)
    return

if __name__ == "__main__":
    main()
