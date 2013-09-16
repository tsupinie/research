
import cPickle
from datetime import datetime, timedelta

import matplotlib
matplotlib.use('agg')
import pylab

import numpy as np

def plotMInnov(minnov_data, times, base_time, styles, colors, labels, ylabel, title, file_name, exclude=[]):
    pylab.figure()
    for exp_name, exp_minnov in minnov_data.iteritems():
        if exp_name not in exclude:
            radar_minnov = np.where(np.isnan(exp_minnov['KCYS']) & (times >= 14400), np.append(exp_minnov['KCYS'][1:], exp_minnov['KCYS'][-1]), exp_minnov['KCYS'])

            pylab.plot(times, radar_minnov, linestyle=styles['KCYS'], color=colors[exp_name], label=labels[exp_name])
#       for radar, radar_minnov in exp_minnov.iteritems():
#           radar_minnov = np.where(np.isnan(radar_minnov) & (times >= 14400), np.append(radar_minnov[1:], radar_minnov[-1]), radar_minnov)
#
#           pylab.plot(times, radar_minnov, linestyle=styles[radar], color=colors[exp_name])

    pylab.legend(loc=4)

    pylab.xticks(times[::2], [ (base_time + timedelta(seconds=int(t))).strftime("%H%M") for t in times[::2] ], rotation=30)
    pylab.xlim([times.min(), times.max()])
    pylab.xlabel("Time (UTC)")

    pylab.ylabel(ylabel)

    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def main():
    base_time = datetime(2009, 6, 5, 18, 0, 0)
    times = np.arange(10800, 18300, 300)[:, np.newaxis].repeat(2, axis=1).flatten('C')

    all_minnov_refl = cPickle.load(open("all_minnov_refl.pkl", 'r'))
    all_minnov_vel = cPickle.load(open("all_minnov_vel.pkl", 'r'))

    styles = {'KCYS':'-', 'KFTG':'--', 'KRIW':':'}
    colors = {'3km-control':'k', '3km-7dBZ,5ms':'r', '3km-n0r=8e5':'g', '3km-fixed-radar':'b', '3kmf-r0h=12km':'m', '3kmf-pr0h=16km':'c', '3kmf-r0h=18km':'#999999'}
    labels = {'3km-control':'CTRL', '3km-7dBZ,5ms':'OBS ERR', '3km-n0r=8e5':"N0R 8e5", '3km-fixed-radar':"CTRL", '3kmf-r0h=12km':"R0H 12", '3kmf-pr0h=16km':"PERT R0H 16", '3kmf-r0h=18km':"R0H 18"}

    plotMInnov(all_minnov_refl, times, base_time, styles, colors, labels, "Mean Innovation (dBZ)", r"Mean Innovation for $Z$ from KCYS (Old Radar Data)", "all_minnov_refl_old_data.png", exclude=['3km-fixed-radar', '3kmf-r0h=12km', '3kmf-pr0h=16km', '3kmf-r0h=18km'])
    plotMInnov(all_minnov_vel, times, base_time, styles, colors, labels, r"Mean Innovation (m s$^{-1}$)", r"Mean Innovation for $v_r$ from KCYS (Old Radar Data)", "all_minnov_vel_old_data.png", exclude=['3km-fixed-radar', '3kmf-r0h=12km', '3kmf-pr0h=16km', '3kmf-r0h=18km'])

    plotMInnov(all_minnov_refl, times, base_time, styles, colors, labels, "Mean Innovation (dBZ)", r"Mean Innovation for $Z$ from KCYS", "all_minnov_refl.png", exclude=['3km-control', '3km-7dBZ,5ms', '3km-n0r=8e5'])
    plotMInnov(all_minnov_vel, times, base_time, styles, colors, labels, r"Mean Innovation (m s$^{-1}$)", r"Mean Innovation for $v_r$ from KCYS", "all_minnov_vel.png", exclude=['3km-control', '3km-7dBZ,5ms', '3km-n0r=8e5'])
    return

if __name__ == "__main__":
    main()
