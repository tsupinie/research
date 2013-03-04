
import numpy as np
from scipy.integrate import odeint

import matplotlib
matplotlib.use('agg')
import pylab
import matplotlib.transforms as transforms
from matplotlib.scale import LogScale

import cPickle

_stlp_data_transform, _stlp_xlabel_transform, _stlp_ylabel_transform = 0, 0, 0
_p_min = 100
_p_max = 1050
_p_step = 100

_T_plot_min = -110
_T_min = -30
_T_max = 50
_T_step = 10
 
_th_min = -30
_th_max = 200
_th_step = 10   

def saturationVaporPressure(temperature):
    L = 2.5e6
    R_v = 461.5
    return 611 * np.exp(L / R_v * (1 / 273.15 - 1 / temperature))

def mixingRatio(vapor_pressure, pressure):
    epsilon = 0.622
    return (epsilon * vapor_pressure) / (pressure - vapor_pressure)

def moistAdiabaticLapseRate(temperature, pressure):
    R_d = 278.
    R_v = 461.5
    c_p = 1003.5
    L = 2.5e6

    sat_mix_rat = mixingRatio(saturationVaporPressure(temperature), pressure)
    moist_term1 = (sat_mix_rat * L) / (R_d * temperature) 
    moist_term2 = (sat_mix_rat * L ** 2) / (R_v * c_p * temperature ** 2)

    return ((1 + moist_term1) * temperature * R_d) / ((1 + moist_term2) * pressure * c_p)    

def pseudoadiabaticLapseRate(temperature, pressure):
    R_d = 278.
    epsilon = 0.622
    c_p = 1003.5
    L = 2.5e6

    sat_mix_rat = mixingRatio(saturationVaporPressure(temperature), pressure)
    moist_term1 = (sat_mix_rat * L) / (R_d * temperature) 
    moist_term2 = (sat_mix_rat * L ** 2 * (epsilon + sat_mix_rat)) / (R_d * c_p * temperature ** 2)
      
    return ((1 + sat_mix_rat) * (1 + moist_term1) * temperature * R_d) / ((1 + sat_mix_rat + moist_term2) * pressure * c_p)    

def _buildTransform(current_axes):
    global _stlp_data_transform, _stlp_xlabel_transform, _stlp_ylabel_transform

    current_figure = current_axes.figure

    current_axes.axes.get_xaxis().set_visible(False)
    current_axes.axes.get_yaxis().set_visible(False)
#   pylab.box(False)

    data_figure_trans = current_axes.transData + current_figure.transFigure.inverted()

    pylab.xlim((_T_min, _T_max))
    pylab.ylim((_p_min, _p_max))

    identity_matrix = np.zeros((3, 3))
    for idx in range(3): identity_matrix[idx, idx] = 1

    # Create the affine matrix for the skew transform.  This only works in data coordinates.  We'll fix that later ...
    skew_matrix = np.copy(identity_matrix)
    skew_matrix[0, 1] = np.tan(45 * np.pi / 180)
    skew_transform = transforms.Affine2D(skew_matrix)

    # Create the logarithmic transform in the y.
    log_p_transform = transforms.blended_transform_factory(transforms.Affine2D(), LogScale(current_axes.yaxis, basey=10).get_transform())

    # The log transform shrinks everything to log(p) space, so define a scale factor to blow it back up to a reasonable size.
    p_bnd_trans = log_p_transform.transform(np.array([[0, _p_min], [0, _p_max]]))[:, 1]
    scale_factor = (_p_max - _p_min) / (p_bnd_trans[1] - p_bnd_trans[0])

    # Define the affine transform for the flip and another for the scale back to reasonable coordinates after the log transform.
    flip_transform = transforms.Affine2D.identity().scale(1, -1)
    preskew_scale_transform = transforms.Affine2D().translate(0, p_bnd_trans[1]).scale(1, scale_factor).translate(0, _p_min)
    postskew_move_transform = transforms.Affine2D().translate(0, _p_min)

    # Define a transform that basically does everything but the skew so we can figure out where the 1000 mb level is and skew around that line.
    prelim_data_transform = log_p_transform + flip_transform + preskew_scale_transform + data_figure_trans
    marker = prelim_data_transform.transform(np.array([[_T_min, 1000]]))[0, 1]

    # Add a translation to that marker point into the data-figure transform matrix.
    data_figure_trans += transforms.Affine2D().translate(0, -marker)

    # Define our skew transform in figure coordinates.
    figure_skew_transform = data_figure_trans + skew_transform + data_figure_trans.inverted()

    # Create our skew-T log-p transform matrix.  It does the log-p transform first, then the flip, then the scale, then the skew.
    _stlp_data_transform = log_p_transform + flip_transform + preskew_scale_transform + figure_skew_transform + current_axes.transData

    # Create a blended transform where the y axis is the log-p, but the x axis is the axes transform (for adding pressure labels and wind barbs).
    _stlp_xlabel_transform = transforms.blended_transform_factory(_stlp_data_transform, current_axes.transAxes)
    _stlp_ylabel_transform = transforms.blended_transform_factory(current_axes.transAxes, _stlp_data_transform)
    return

def plotSkewTBackground(current_axes):
    current_axes.set_position([0, 0, 1, 1])

    current_figure = current_axes.figure
    current_figure.set_dpi(100)

    if type(_stlp_data_transform) == int:
        _buildTransform(current_axes)

    # Draw the isobars.
    for p_line in range(_p_min, _p_max + _p_step, _p_step):
        pylab.plot([_T_plot_min, _T_max], [p_line, p_line], color='k', lw=0.75, transform=_stlp_data_transform)
        pylab.text(0, p_line, "%d" % p_line, va='center', transform=_stlp_ylabel_transform)

    # Draw the isotherms.
    for T_line in range(_T_plot_min, _T_max + _T_step, _T_step):
        if T_line < 0:
            weight = 0.75
            color = 'c'
        elif T_line > 0: 
            color = '#880000'
            weight = 0.75
        else: 
            weight = 1
            color = 'b'

        pylab.plot([T_line, T_line], [_p_min, _p_max], color=color, lw=weight, transform=_stlp_data_transform)
        pylab.text(T_line, 1000, "%d" % T_line, ha='center', transform=_stlp_data_transform)

    # Draw the dry adiabats.
    for th_line in range(_th_min, _th_max + _th_step, _th_step):
        p_th = np.arange(_p_min, _p_max + 10., 10.)
        t_th = (th_line + 273.15) * (p_th / 1000.) ** (2. / 7.) - 273.15
        pylab.plot(t_th, p_th, color='#802a2a', lw=0.75, transform=_stlp_data_transform)

    # Draw the mixing ratio lines.
    for w_line in [ 1, 2, 3, 5, 8, 12, 16, 20 ]:
        p_w = np.arange(600., _p_max + 10., 10.)
        e_w = (p_w * 100 * w_line / 1000.) / (w_line / 1000. + 0.622)
        td_w = 1 / (1/273.0 - 461.5 / 2.5e6 * np.log(e_w / 611.)) - 273.15
        pylab.plot(td_w, p_w, color='#7fff00', lw=0.75, linestyle='--', transform=_stlp_data_transform)
        pylab.text(td_w[0], p_w[0], "%d" % w_line, color="#7fff00", ha='center', transform=_stlp_data_transform)

    # Draw the moist adiabats.
    for the_line in xrange(-20, 55, 5):
        p_the_above = np.arange(_p_min, 1010., 10.)
        p_the_below = np.arange(1000., _p_max + 10., 10.)

        t_the_above = odeint(moistAdiabaticLapseRate, the_line + 273.15, 100 * p_the_above[::-1])[::-1,0] - 273.15
        t_the_below = odeint(moistAdiabaticLapseRate, the_line + 273.15, 100 * p_the_below)[:,0] - 273.15

        p_the = np.concatenate((p_the_above, p_the_below[1:]))
        t_the = np.concatenate((t_the_above, t_the_below[1:]))

        pylab.plot(t_the, p_the, color='#2e8b57', lw=0.75, linestyle=':', transform=_stlp_data_transform)

#   pylab.plot([0.95, 0.95], [p_min, p_max], color='k', lw=0.5, transform=stlp_ylabel_transform)

#   t_barb_locs, p_barb_locs = np.swapaxes(stlp_ylabel_transform.transform(np.dstack((0.95 * np.ones(good_idxs.shape), p_obs_snd[good_idxs]))[0]), 0, 1)

#   print t_barb_locs
#   print p_barb_locs

#   pylab.barbs(t_barb_locs, p_barb_locs, u_obs_snd[good_idxs], v_obs_snd[good_idxs])

    pylab.xlim((_T_min, _T_max))
    pylab.ylim((_p_min, _p_max))
    return

def plotProfile(t_snd, p_snd, **kwargs):
    kwargs['transform'] = _stlp_data_transform
    plot_break = 90 # For some reason, one sounding gets cut off if you try to plot more than 100 values at the same time ... ?????    

    pylab.plot(t_snd[:plot_break], p_snd[:plot_break], **kwargs)
    pylab.plot(t_snd[(plot_break - 1):], p_snd[(plot_break - 1):], **kwargs)

    pylab.xlim((_T_min, _T_max))
    pylab.ylim((_p_min, _p_max))
    return

def plotWinds(u_snd, v_snd, p_snd, **kwargs):
    t_base = 45
    stride = 50
    thin_data = (slice(None, None, stride))
    pylab.axvline(x=t_base, color='k', lw=0.5)

    points = np.vstack((t_base * np.ones(p_snd.shape), p_snd)).T

    x, base_y = _stlp_data_transform.transform(np.array([[t_base, 1000]]))[0]
    trans_points = _stlp_data_transform.transform(points)
    trans_points[:, 0] = x

    axes_points = pylab.gca().transData.inverted().transform(trans_points)

    plot_xs, plot_ys = zip(*axes_points)
    pylab.barbs(plot_xs[thin_data], plot_ys[thin_data], u_snd[thin_data], v_snd[thin_data])
    return

def plotSounding(file_name, **kwargs):
    pylab.clf()
    plotSkewTBackground(pylab.gca())

    plotProfile(kwargs['t'], kwargs['p'], color='r', lw=1.5)
    plotProfile(kwargs['td'], kwargs['p'], color='g', lw=1.5)

    pylab.savefig(file_name)
    return

def main():
    soundings = cPickle.load(open('soundings.pkl', 'r'))
    for sounding in soundings:
        u_obs_snd = sounding['u_wind']
        v_obs_snd = sounding['v_wind']
        t_obs_snd = sounding['temperature']
        td_obs_snd = sounding['dewpoint']
        p_obs_snd = sounding['pressure']

        good = (u_obs_snd != 9999.0) & (v_obs_snd != 9999.0) & (t_obs_snd != 999.0) & (td_obs_snd != 999.0) & (p_obs_snd != 9999.0)
        good_idxs = np.where(good)[0]

        name = sounding['release_site'][:4]
        name = name.replace('/', '_')
        plotSounding("proximity.%s.png" % name, p=p_obs_snd[good_idxs], t=t_obs_snd[good_idxs], td=td_obs_snd[good_idxs], u=u_obs_snd[good_idxs], v=v_obs_snd[good_idxs])

    return
if __name__ == "__main__":
    main()
