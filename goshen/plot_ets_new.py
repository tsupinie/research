
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

import cPickle
from datetime import datetime, timedelta

def plotETS(ets_data, times, base_time, styles, colors, labels, title, file_name, exclude=[]):
    pylab.figure()
    for exp_name in sorted(ets_data.keys()):
        exp_ets = ets_data[exp_name]
        if exp_name not in exclude:
            pylab.plot(times, exp_ets['KCYS'], linestyle=styles['KCYS'], color=colors[exp_name], label=labels[exp_name])
#       for radar, radar_ets in exp_ets.iteritems():
#           pylab.plot(times, radar_ets, linestyle=styles[radar], color=colors[exp_name])

    pylab.legend()

    pylab.xticks(times, [ (base_time + timedelta(seconds=int(t))).strftime("%H%M") for t in times ], rotation=30)
    pylab.xlim([times.min(), times.max()])
    pylab.xlabel("Time (UTC)")

    pylab.ylabel("ETS")
    pylab.ylim([0., 0.9])

    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def main():
    base_time = datetime(2009, 6, 5, 18, 0, 0)
    times = np.arange(14400, 18300, 300)
    refl_threshold = 40
    vel_threshold = 30

    all_ets_refl = cPickle.load(open("all_ets_new_1km_%02ddBZ_bnd.pkl" % refl_threshold, 'r'))
    all_ets_vel = cPickle.load(open("all_ets_new_1km_%02dms_bnd.pkl" % vel_threshold, 'r'))

    styles_3km = {'KCYS':'-', 'KFTG':'--', 'KRIW':':'}
    colors_3km = {'3km-control':'k', '3km-7dBZ,5ms':'r', '3km-n0r=8e5':'g', '3km-fixed-radar':'b', '3kmf-r0h=12km':'m', '3kmf-pr0h=16km':'c', '3kmf-r0h=18km':'#999999'}
    labels_3km = {'3km-control':'CTRL', '3km-7dBZ,5ms':'OBS ERR', '3km-n0r=8e5':"N0R 8e5", '3km-fixed-radar':"CTRL", '3kmf-r0h=12km':"R0H 12", '3kmf-pr0h=16km':"PERT R0H 16", '3kmf-r0h=18km':"R0H 18"}

    styles_1km = {'KCYS':'-', 'KFTG':'--', '05XP':':'}
    colors_1km = {'1kmf-zupdtpt':'k', '1kmf-z-no-05XP':'r', '1kmf-z-no-mm-05XP':'g', '1kmf-z-no-mm':'b', '1kmf-z-no-v2':'m', '1kmf-z-no-snd':'c' }
    labels_1km = {'1kmf-zupdtpt':'CTRL', '1kmf-z-no-05XP':'NO MWR05XP', '1kmf-z-no-mm-05XP':"NO MM OR MWR05XP", '1kmf-z-no-mm':"NO MM", '1kmf-z-no-v2':"NO V2 OBS", '1kmf-z-no-snd':"NO SND"}

    styles = styles_1km
    colors = colors_1km
    labels = labels_1km

#   plotETS(all_ets_refl, times, base_time, styles, colors, labels, r"ETS for $Z$ > %02d dBZ from KCYS (Old radar data)" % refl_threshold, "all_ets_new_1km_%02ddBZ_old_data.png" % refl_threshold)
#   plotETS(all_ets_vel, times, base_time, styles, colors, labels, r"ETS for |$v_r$| > %02d m s$^{-1}$ from KCYS (Old radar data)" % vel_threshold, "all_ets_new_1km_%02dms_old_data.png" % vel_threshold)

    plotETS(all_ets_refl, times, base_time, styles, colors, labels, r"ETS for $Z$ > %02d dBZ from KCYS" % refl_threshold, "all_ets_new_1km_%02ddBZ_bnd.png" % refl_threshold)
    plotETS(all_ets_vel, times, base_time, styles, colors, labels, r"ETS for |$v_r$| > %02d m s$^{-1}$ from KCYS" % vel_threshold, "all_ets_new_1km_%02dms_bnd.png" % vel_threshold)
    return

if __name__ == "__main__":
    main()
