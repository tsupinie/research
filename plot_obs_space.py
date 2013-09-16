
import numpy as np

from assemble_rmsd import _loadRMSDFile, _loadSpreadFile, _loadconsistencyparameterfile
from temporal import goshen_1km_temporal
from util import publicationFigure

from collections import OrderedDict

import matplotlib
#matplotlib.use('agg')
import pylab

params = {
    'innov': ('rmsd',   {'anal':'an',  'fcst':'bg'},  _loadRMSDFile),
    'spread':('spread', {'anal':'ana', 'fcst':'fcs'}, _loadSpreadFile),
    'cons':  ('innov',  {'anal':'',    'fcst':''},    _loadconsistencyparameterfile),
}

def sawtooth(fcst, anal):
    return np.vstack((fcst, anal)).flatten('F')

def loadRadar(path, radar, temp, quants=['vr', 'Z'], states=['fcst', 'anal']):
    radar_params = {}
    for param, (file_base, state_map, load_fn) in params.iteritems():
        sawtooth_data = []
        for state in states:
            radar_quants = dict( (q, []) for q in quants )

            for time in temp:
                file_name = "%s/%s%s%s%06d" % (path, radar, file_base, state_map[state], time)
                for q, dat in zip(quants, load_fn(file_name)):
                    radar_quants[q].append(dat)

            sawtooth_data.append(radar_quants)
        radar_params[param] = dict( (q, sawtooth(*[ s[q] for s in sawtooth_data ])) for q in quants )
    return radar_params

def plotRadar(rad_data, temp, exp_names, quant, file_name):
    # rad_data: experiment, parameter, observed quantity
    colors = dict(zip(exp_names.keys(), ['k', 'r', 'g', 'b', 'c', 'm']))

    def consistencySubplot(multiplier=1.0, layout=(-1, -1)):
        cons_rat_data = dict( (e, rad_data[e]['cons'][quant]) for e in rad_data.iterkeys() )

        times = temp.getTimes()

        for exp, exp_name in exp_names.iteritems():
            pylab.plot(sawtooth(times, times), cons_rat_data[exp], color=colors[exp], label=exp_name)

        pylab.axhline(y=1, color='k', linestyle=':')
        pylab.axvline(x=14400, color='k', linestyle=':')

        pylab.xlabel("Time (UTC)", size='large')
        pylab.ylabel("Consistency Ratio", size='large')

        pylab.xlim(times[0], times[-1])

        pylab.xticks(times[::2], temp.getStrings("%H%M", aslist=True)[::2], rotation=30, size='x-large')
        pylab.yticks(size='x-large')
        pylab.legend()
        return

    def rmsdSpreadSubplot(multiplier=1.0, layout=(-1, -1)):
        rmsd_data   = dict( (e, rad_data[e]['innov'][quant])  for e in rad_data.iterkeys() )
        spread_data = dict( (e, rad_data[e]['spread'][quant]) for e in rad_data.iterkeys() )

        times = temp.getTimes()
        n_t = len(times)

        for exp, exp_name in exp_names.iteritems():
            pylab.plot(sawtooth(times, times)[:(n_t + 1)], rmsd_data[exp][:(n_t + 1)], color=colors[exp], linestyle='-')
            pylab.plot(times[(n_t / 2):], rmsd_data[exp][n_t::2], color=colors[exp], linestyle='-')
 
        for exp, exp_name in exp_names.iteritems():
            pylab.plot(sawtooth(times, times)[:(n_t + 1)], spread_data[exp][:(n_t + 1)], color=colors[exp], linestyle='--')
            pylab.plot(times[(n_t / 2):], spread_data[exp][n_t::2], color=colors[exp], linestyle='--')

        ylim = pylab.ylim()
        pylab.plot(times, -1 * np.ones((len(times),)), color='#999999', linestyle='-', label="RMS Innovation")
        pylab.plot(times, -1 * np.ones((len(times),)), color='#999999', linestyle='--', label="Spread")

        pylab.axhline(y=7, color='k', linestyle=':')
        pylab.axvline(x=14400, color='k', linestyle=':')

        pylab.ylabel("RMS Innovation/Spread (dBZ)", size='large')

        pylab.xlim(times[0], times[-1])
        pylab.ylim(ylim)

        pylab.legend(loc=4)

        pylab.xticks(times[::2], [ "" for t in times[::2] ])
        pylab.yticks(size='x-large')
        return

    pylab.figure(figsize=(8, 9.5))
    pylab.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1, hspace=0.05)
    publicationFigure([rmsdSpreadSubplot, consistencySubplot], (2, 1), corner='ul')

    pylab.savefig(file_name)
    pylab.close()
    return

def main():
    base_path = "/caps2/tsupinie/"
    radars = [ 'KCYS' ] #, 'KFTG, '05XP' ]
    experiments = OrderedDict([('1kmf-sndr0h=25km', 'CTRL'), ('1kmf-zs25-no-05XP', 'NO_MWR'), ('1kmf-z-no-snd', 'NO_SND'), ('1kmf-zs25-no-mm', 'NO_MM'), ('1kmf-zs25-no-mm-05XP', 'NO_MWR_MM'), ('1kmf-z-no-v2', 'NO_V2')])
    temp = goshen_1km_temporal()

    for rad in radars:
        radar_params = {}

        for exp, name in experiments.iteritems():
            radar_params[exp] = loadRadar("%s/%s" % (base_path, exp), rad, temp)

        plotRadar(radar_params, temp, experiments, 'Z', "obs_space_%s.png" % rad)
    return

if __name__ == "__main__":
    main()
