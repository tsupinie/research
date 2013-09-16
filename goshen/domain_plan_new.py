
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
from matplotlib import transforms

from grid import goshen_1km_grid, goshen_3km_grid
from temporal import goshen_1km_temporal, goshen_3km_temporal
from radarobsfile import RadarObsFile
from util import publicationFigure

import glob
import string

#def publicationFigure(subfigures, layout, corner='ul'):
#    rows, columns = layout
#
#    bbox = {'ec':'k', 'fc':'w', 'pad':10, 'lw':2}
#    loc = {'ul':(0.0, 1.0)}
#    align = {'ul':('left', 'top')}
#    for idx, sf in enumerate(subfigures):
#        pylab.subplot(rows, columns, idx + 1)
#        sf()
#
#        text_x, text_y = loc[corner]
#        h_align, v_align = align[corner]
#        pylab.text(text_x, text_y, "(%s)" % string.ascii_lowercase[idx], transform=pylab.gca().transAxes, bbox=bbox, ha=h_align, va=v_align, fontsize=32, fontweight='bold')
#    return

def main():
    radar_path = "/data6/tsupinie/goshen/qc/1km/"
    nom_2_real = {'2100':'2058', '2200':'2158', '2300':'2257'}

    radar_locations = {
        'KCYS':((41.15194, -104.80611), 230),
        'KFTG':((39.78667, -104.54583), 230),
        'KRIW':((43.06611, -108.47722), 230),
        'MWR-05XP':((41.56150, -104.298996), 40),
    }

    boxes = (
        (slice(65, 105),  slice(115, 155)),
        (slice(105, 145), slice(115, 155)),
        (slice(140, 180), slice(110, 150)),
    )

    refl = tuple(RadarObsFile(f)['Z'][0] for f in sorted(glob.glob("%s/KCYS.20090605.2[123]00*" % radar_path)))

    buffer_1km = 0
    grid_1km_plot = goshen_1km_grid(buffer=buffer_1km)
    grid_1km = goshen_1km_grid()
    grid_1km_lats, grid_1km_lons = grid_1km.getBoundaryCoords()

    grid_3km_plot = goshen_3km_grid()#buffer=40)
    grid_3km = goshen_3km_grid()
    grid_3km_lats, grid_3km_lons = grid_3km.getBoundaryCoords()

    temp = goshen_1km_temporal(step=3600)
    temp_1km = goshen_1km_temporal()
    temp_3km = goshen_3km_temporal()

    def sf_3km_domain(multiplier=1.0, layout=(-1, -1)):
        pylab.plot(*grid_3km_plot(grid_3km_lons, grid_3km_lats), color='k', lw=multiplier)
        pylab.plot(*grid_3km_plot(grid_1km_lons, grid_1km_lats), color='k', lw=2 * multiplier)

        width, height = grid_3km_plot.getWidthHeight()
        pylab.xlim(0, width)
        pylab.ylim(0, height)

        for radar_id, ((radar_lat, radar_lon), radar_range) in radar_locations.iteritems():
            radar_x, radar_y = grid_3km_plot(radar_lon, radar_lat)

            if radar_id == "MWR-05XP":
                offset_sign = 1
                h_align, v_align = 'left', 'bottom'
            else:
                offset_sign = -1
                h_align, v_align = 'right', 'top'

            ax_trans = pylab.gca().transData + pylab.gca().transAxes.inverted() 
            label_x, label_y = ax_trans.transform(np.array([[radar_x, radar_y]]))[0] + offset_sign * 0.0125

            if radar_id not in [ 'MWR-05XP' ]:
                pylab.plot(radar_x, radar_y, 'k^', ms=6 * multiplier)
                pylab.text(label_x, label_y, radar_id, color='k', clip_box=pylab.gca().bbox, size=12*multiplier, clip_on=True, ha=h_align, va=v_align, transform=pylab.gca().transAxes)
                pylab.gca().add_patch(Circle((radar_x, radar_y), 1000 * radar_range, fc='none', lw=multiplier, zorder=10))

        grid_3km_plot.drawPolitical(color='#999999', lw=1.5 * multiplier)

    def sf_1km_domain(multiplier=1.0, layout=(-1, -1)):
        pylab.plot(*grid_1km_plot(grid_1km_lons, grid_1km_lats), color='k', lw=multiplier)

        width, height = grid_1km_plot.getWidthHeight()
        x, y = grid_1km_plot.getXY()

        for r, b, t in zip(refl, boxes, temp.getStrings("%H%M")):
            bx, by = b

            plot_r = np.zeros(r.shape, dtype=r.dtype)
            plot_r[b[::-1]] = r[b[::-1]]

            if buffer_1km == 0:
                buffer_slc = slice(None, None)
            else:
                buffer_slc = slice(buffer_1km, buffer_1km)

            pylab.contour(x[buffer_slc, buffer_slc], y[buffer_slc, buffer_slc], plot_r, levels=[ 40. ], colors='k', linewidths=multiplier, zorder=10)
            pylab.text((bx.stop + buffer_1km) * 1000, (by.stop + buffer_1km) * 1000, nom_2_real[t], color='k', size=12*multiplier, ha='right', va='top')

        pylab.xlim(0, width)
        pylab.ylim(0, height)

        for radar_id, ((radar_lat, radar_lon), radar_range) in radar_locations.iteritems():
            radar_x, radar_y = grid_1km_plot(radar_lon, radar_lat)

            offset_sign = -1

            ax_trans = pylab.gca().transData + pylab.gca().transAxes.inverted()
            label_x, label_y = ax_trans.transform(np.array([[radar_x, radar_y]]))[0] + offset_sign * 0.0125

            if radar_id not in [ 'KRIW' ]:
                pylab.plot(radar_x, radar_y, 'k^', ms=6 * multiplier)
                pylab.text(label_x, label_y, radar_id, color='k', clip_box=pylab.gca().bbox, clip_on=True, size=12*multiplier, ha='right', va='top', transform=pylab.gca().transAxes)
                pylab.gca().add_patch(Circle((radar_x, radar_y), 1000 * radar_range, fc='none', lw=multiplier, zorder=10))

        grid_1km_plot.drawPolitical(color='#999999', lw=1.5 * multiplier)

    pylab.figure(figsize=(16, 12))
    pylab.subplots_adjust(left=0.05, bottom=0.375, right=0.95, top=0.975, hspace=0.1, wspace=0.1)
    publicationFigure([ sf_3km_domain, sf_1km_domain ], (1, 2))

    tl_ax = pylab.axes((0.05, 0.025, 0.9, 0.325))
    tl_ax.xaxis.set_ticks([])
    tl_ax.yaxis.set_visible(False)

    center_3km = 0.65
    center_1km = 0.35
    rad_spacing = 0.05

    tl_ax.arrow(0, center_3km, 18300, 0, head_width=0.025, head_length=300, fc='k')
    pylab.text(-300, center_3km, "Outer Domain", ha='right', va='center', size=18)
    for t_ens, t_str in zip(temp_3km, temp_3km.getStrings("%H%M")):
        if t_ens % 3600:
            tick_size = 0.0125
        else:
            tick_size = 0.025
            pylab.text(t_ens, 0.9, "%s UTC" % t_str, size=22, ha='center', va='center')
        tl_ax.add_line(Line2D([t_ens, t_ens], [center_3km - tick_size, center_3km + tick_size], color='k', lw=2))

    for idx, (t_st, t_en, rad, c) in enumerate([[1800, 14400, 'KRIW', 'r'], [5400, 14400, 'KCYS', 'g'], [5400, 14400, 'KFTG', 'b']]):
        tl_ax.add_line(Line2D([t_st, t_en], [ center_3km + (idx + 1) * rad_spacing ] * 2, color=c, lw=5))
        pylab.text(t_st - 300, center_3km + (idx + 1) * rad_spacing, rad, ha='right', va='center', color=c, size=18)

    tl_ax.arrow(10800, center_1km, 7500, 0, head_width=0.025, head_length=300, fc='k')
    pylab.text(-300, center_1km, "Inner Domain", ha='right', va='center', size=18)
    for t_ens in temp_1km:
        if t_ens % 3600:
            tick_size = 0.0125
        else:
            tick_size = 0.025
        tl_ax.add_line(Line2D([t_ens, t_ens], [center_1km - tick_size, center_1km + tick_size], color='k', lw=2))

    for idx, (t_st, t_en, rad, c) in enumerate([[11100, 14400, 'KCYS', 'g'], [11700, 14400, 'KFTG', 'b'], [13500, 14400, 'MWR-05XP', 'k']]):
        tl_ax.add_line(Line2D([t_st, t_en], [ center_1km - (idx + 1) * rad_spacing ] * 2, color=c, lw=5))
        pylab.text(t_st - 300, center_1km - (idx + 1) * rad_spacing, rad, ha='right', va='center', color=c, size=18)

    tl_ax.arrow(10800, center_3km - rad_spacing, 0, -(center_3km - rad_spacing) + (center_1km + rad_spacing), head_width=150, head_length=0.025, length_includes_head=True, fc='k')
    tl_ax.add_line(Line2D([13920, 16230], [0.5, 0.5], color='m', lw=5))
    pylab.text(13620, 0.5, "Tornado", ha='right', va='center', size=18, color='m')

    base_diag = np.hypot(6, 12)
    size_x, size_y = pylab.gcf().get_size_inches()
    fig_diag = np.hypot(size_x, size_y)

    multiplier = fig_diag / base_diag

    corner = 'ul'
    bbox = {'ec':'k', 'fc':'w', 'pad':10, 'lw':multiplier}
    loc = {'ul':(0.0, 1.0), 'ur':(1.0, 1.0)}
    align = {'ul':('left', 'top'), 'ur':('right', 'top')}
    pad_signs = {'l':1, 'r':-1, 'u':-1}

    offset = transforms.ScaledTranslation(pad_signs[corner[1]] * ((bbox['pad'] + bbox['lw']) / (2 * 72.)), pad_signs[corner[0]] * (bbox['pad'] + bbox['lw']) / (2 * 72.), pylab.gcf().dpi_scale_trans)
    text_transform = pylab.gca().transAxes + offset
    text_x, text_y = loc[corner]
    h_align, v_align = align[corner]
    pylab.text(text_x, text_y, "(c)", transform=text_transform, bbox=bbox, ha=h_align, va=v_align, fontsize=16 * multiplier, fontweight='bold', zorder=1000)

    pylab.xlim(-3300, 19200)
    pylab.ylim(0, 1)

    pylab.savefig("domain_plan.png")
    pylab.close()

    return

if __name__ == "__main__":
    main()
