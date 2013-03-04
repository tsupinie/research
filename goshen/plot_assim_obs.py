
from util import goshen_1km_proj, goshen_1km_gs, goshen_3km_proj, goshen_3km_gs, setupMapProjection, drawPolitical

import cPickle

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon, Circle
from matplotlib.gridspec import GridSpec

def plotSoundingObservationsComposite(obs, map, scale_len, dimensions, title, file_name, radars=None):
    figure = pylab.figure(figsize=(10,8))
    x_size, y_size = dimensions

    grid_spec = GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 3])

    colors = [ 'r', 'g', 'b', 'c', 'm', '#660099', '#006666', '#99ff00', '#ff9900', '#666666' ]
    names = [ "Sounding (NSSL)", "Sounding (NCAR)", "Sounding (NCAR)", "Sounding (NCAR)", "Sounding (NSSL)", "Sounding (NSSL)" ]
    ob_ids = np.unique1d(obs['id'])

    bottom = 0.075
    top = 0.95
    stretch_fac = (1 - top + bottom) / 2. - 0.011

    pylab.subplots_adjust(left=0.1 + stretch_fac, bottom=bottom, right=0.9 - stretch_fac, top=top, hspace=0.075, wspace=0.075)

    ax_xy = figure.add_subplot(grid_spec[2])
    ax_xz = figure.add_subplot(grid_spec[0], sharex=ax_xy)
    ax_yz = figure.add_subplot(grid_spec[3], sharey=ax_xy)

    def labelMe(do_label, label):
        if do_label:
            return not do_label, label
        else:
            return do_label, None

    color_num = 0
    prof_label = True
    sfc_label = True
    asos_label = True
    ttu_label = True

    for ob_id in ob_ids:
        plot_vertical = False
        ob_idxs = np.where(obs['id'] == ob_id)[0]
        these_obs = obs[ob_idxs]
        ob_xs, ob_ys = map(these_obs['longitude'], these_obs['latitude'])

        if np.any(ob_xs >= 0) and np.any(ob_xs <= x_size) and np.any(ob_ys >= 0) and np.any(ob_ys <= y_size):
            if these_obs['obtype'][0] == "PROF":
                prof_label, label = labelMe(prof_label, "Profiler")
                color = "#ff9900"
                marker = 'o'
                ms = 3
                plot_vertical = True
            elif these_obs['obtype'][0] == "SNDG":
                color = colors[color_num]
                label = "Sounding"# % ob_id
                marker = 'o'
                ms = 3
                plot_vertical = True
                color_num += 1               
            elif ob_id[0] == "P":
                color = "#999999"
                sfc_label, label = labelMe(sfc_label, "Sounding")
                marker = 'o'
                ms = 3
            elif ob_id[0] == "K":
                asos_label, label = labelMe(asos_label, "ASOS")
                color = 'k'
                marker = '*'
                ms = 5
            elif ob_id[0] == "1" or ob_id[0] == "2":
                ttu_label, label = labelMe(ttu_label, "TTU Sticknet")
                color = "#003300"
                marker = 'o'
                ms = 3

            ax_xy.plot(ob_xs, ob_ys, marker, mfc=color, mec=color, ms=ms, label=label)[0]
            if plot_vertical:
                ax_xz.plot(ob_xs, these_obs['elevation'] / 1000., marker, mfc=color, mec=color, ms=ms)
                ax_yz.plot(these_obs['elevation'] / 1000., ob_ys, marker, mfc=color, mec=color, ms=ms)

    dummy, zlim = ax_xz.get_ylim()

    for ob_id in ob_ids:
        ob_idxs = np.where(obs['id'] == ob_id)[0]
        these_obs = obs[ob_idxs]
        ob_xs, ob_ys = map(these_obs['longitude'], these_obs['latitude'])

        plot_vertical = False

        if np.any(ob_xs >= 0) and np.any(ob_xs <= x_size) and np.any(ob_ys >= 0) and np.any(ob_ys <= y_size):
            if these_obs['obtype'][0] == "PROF":
                color = "#ff9900"
                plot_vertical = True
            elif these_obs['obtype'][0] == "SNDG":
                color = "#999999"
                plot_vertical = True

            if plot_vertical:
                ax_xy.plot(ob_xs[0], ob_ys[0], 'x', color=color, ms=6)
                ax_xz.plot([ob_xs[0], ob_xs[0]], [0, zlim], '--', color=color, lw=0.5)
                ax_yz.plot([0, zlim], [ob_ys[0], ob_ys[0]], '--', color=color, lw=0.5)

    if radars:
        rad_label = True
        for radar in radars:
            rad_label, label = labelMe(rad_label, "Radar")
            (lat, lon), range = radar
            radar_x, radar_y = map(lon, lat)
            ax_xy.plot(radar_x, radar_y, 'b<', ms=4, label=label)
            ax_xy.add_patch(Circle((radar_x, radar_y), range, fc='none', ec='k'))

    ax_xz.set_ylim(0, zlim)
    ax_yz.set_xlim(0, zlim)

    pylab.sca(ax_xy)
    drawPolitical(map, scale_len=scale_len)

    ax_xz.set_xlabel("x")
    ax_xz.set_ylabel("z (km)")
    ax_yz.set_xlabel("z (km)")
    ax_yz.set_ylabel("y", rotation=-90)

    line_objects = [ l for l in sorted(ax_xy.lines, key=lambda x: x.get_label()) if l.get_label()[0] != "_" ]
    ax_xy.legend(line_objects, [ l.get_label() for l in line_objects], loc=1, numpoints=1, ncol=2, prop={'size':'medium'}, bbox_to_anchor=(1.28, 1.07, 0.33, 0.33), handletextpad=0, columnspacing=0)

    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def plotObservationsComposite(obs, map, scale_len, title, file_name):
    pylab.figure(figsize=(10,8))
    pylab.axes((0, 0, 1, 0.95))
    colors = [ 'r', 'g', 'b', 'c', 'm', '#660099', '#ff9900', '#006666' ]

    ob_ids = np.unique1d(obs['id'])
    ttu_label = False
    asos_label = False
    sndg_label = False

    for ob_id in ob_ids:
        ob_idxs = np.where(obs['id'] == ob_id)[0]
        these_obs = obs[ob_idxs]
        ob_xs, ob_ys = map(these_obs['longitude'], these_obs['latitude'])

        if ob_id[0] == "P":
            ob_num = int(ob_id[1]) - 1
            pylab.plot(ob_xs, ob_ys, 'o', mfc=colors[ob_num], mec=colors[ob_num], ms=3, label="NSSL MM (%s)" % ob_id)
        elif ob_id[0] == "K":
            if not asos_label:
                label = "ASOS"
                asos_label = True
            else:
                label = None

            pylab.plot(ob_xs[0], ob_ys[0], 'k*', ms=5, label=label)
        elif ob_id[0] == "1" or ob_id[0] == "2":
            if not ttu_label:
                label = "TTU Sticknet"
                ttu_label = True
            else:
                label = None

            pylab.plot(ob_xs[0], ob_ys[0], 'o', mfc="#999999", mec="#999999", ms=3, label=label)
        elif these_obs[0]['obtype'] == "SNDG":
            if not sndg_label:
                label = "Sounding"
                sndg_label = True
            else:
                label = None

            pylab.plot(ob_xs[0], ob_ys[0], 'k^', ms=4, label=label)

    drawPolitical(map, scale_len=scale_len)

    line_objects = [ l for l in sorted(pylab.gca().lines, key=lambda x: x.get_label()) if l.get_label()[0] != "_" ]
    pylab.legend(line_objects, [ l.get_label() for l in line_objects], loc=2, numpoints=1, prop={'size':'medium'})
    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def main():
    domain = "3km"

    assim_obs = cPickle.load(open("assim_obs_%s.pkl" % domain, 'r'))
    sounding_obs = cPickle.load(open("sounding_obs_da_%s.pkl" % domain, 'r'))

    radar_location_KCYS = (41.15194, -104.80611)
    radar_location_KFTG = (39.78667, -104.54583)
    radar_location_KRIW = (43.06611, -108.47722)
    radar_location_05XP = (41.56150, -104.298996)

    if domain == "3km":
        time_span = "1830-2200 UTC"
        goshen_proj = goshen_3km_proj
        goshen_gs = goshen_3km_gs
        radars = [ (radar_location_KCYS, 230000), (radar_location_KFTG, 230000), (radar_location_KRIW, 230000) ] 

        profiler_obs = cPickle.load(open("profile_obs_da.pkl", 'r'))
        scale = 75
    elif domain == "1km":
        time_span = "2100-2200 UTC"
        goshen_proj = goshen_1km_proj
        goshen_gs = goshen_1km_gs
        radars = [ (radar_location_KCYS, 230000), (radar_location_KFTG, 230000), (radar_location_05XP, 40000) ] 

        profiler_obs = np.empty((0,), dtype=sounding_obs.dtype)
        scale = 25

    sounding_obs = np.append(sounding_obs, profiler_obs)

    proj = setupMapProjection(goshen_proj, goshen_gs)
    map = Basemap(**proj)

#   assim_obs['sndg'] = sounding_obs

    plotObservationsComposite(np.concatenate(tuple(assim_obs.values())), map, scale, "Surface Observation Composite (%s)" % time_span, "mm_composite_%s.png" % domain)
    plotSoundingObservationsComposite(sounding_obs, map, scale, (goshen_proj['width'], goshen_proj['height']), 
        "Sounding Observation Composite (%s)"  % time_span, "sndg_composite_%s.png" % domain, radars=radars)
    return

if __name__ == "__main__":
    main()
