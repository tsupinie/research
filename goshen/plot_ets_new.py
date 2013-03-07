
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

import cPickle

def plotETS(ets_data, times, styles, colors, title, file_name):
    pylab.figure()
    for exp_name, exp_ets in ets_data.iteritems():
        for radar, radar_ets in exp_ets.iteritems():
            pylab.plot(times, radar_ets, linestyle=styles[radar], color=colors[exp_name])

    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def main():
    times = np.arange(14400, 18300, 300)
    refl_threshold = 20
    vel_threshold = 20

    all_ets_refl = cPickle.load(open("all_ets_new_%02ddBZ.pkl" % refl_threshold, 'r'))
    all_ets_vel = cPickle.load(open("all_ets_new_%02dms.pkl" % vel_threshold, 'r'))

    styles = {'KCYS':'-', 'KFTG':'--', 'KRIW':':'}
    colors = {'3km-control':'k', '3km-7dBZ,5ms':'r', '3km-n0r=8e5':'g', '3km-fixed-radar':'b', '3kmf-r0h=12km':'m', '3kmf-pr0h=16km':'c', '3kmf-r0h=18km':'#999999'}

    plotETS(all_ets_refl, times, styles, colors, r"ETS for $Z$ > %02d dBZ" % refl_threshold, "all_ets_new_%02ddBZ.png" % refl_threshold)
    plotETS(all_ets_vel, times, styles, colors, r"ETS for |$v_r$| > %02d m s$^{-1}$" % vel_threshold, "all_ets_new_%02dms.png" % vel_threshold)

    return

if __name__ == "__main__":
    main()
