
import numpy as np
import matplotlib
matplotlib.use('agg')
import pylab

from datetime import datetime, timedelta
import cPickle

radar_id_1km_goshen = [ "KTLX" ]
experiments_cov_infl = [ "mult=1.15", "mult=1.15,noise=0.5", "mult=1.15,noise=0.75", "adapt=0.80", "adapt=0.80,noise=0.5", "adapt=0.80,noise=0.75", "relax=0.50", "relax=0.50,noise=0.5", "relax=0.50,noise=0.75", "no-infl" ]
#experiments_cov_infl = [ "mult=1.15", "adapt=0.80", "adapt=0.80,noise=0.5", "adapt=0.80,noise=0.75", "relax=0.50", "relax=0.50,noise=0.5", "relax=0.50,noise=0.75", "no-infl" ]
colors_cov_infl = {
    "mult=1.15:KTLX":'r',  "mult=1.15,noise=0.5:KTLX":'r',  "mult=1.15,noise=0.75:KTLX":'r', 
    "adapt=0.80:KTLX":'g', "adapt=0.80,noise=0.5:KTLX":'g', "adapt=0.80,noise=0.75:KTLX":'g', 
    "relax=0.50:KTLX":'b', "relax=0.50,noise=0.5:KTLX":'b', "relax=0.50,noise=0.75:KTLX":'b',
    "no-infl:KTLX":'m'
}
styles_cov_infl = {
    "mult=1.15:KTLX":'-',  "mult=1.15,noise=0.5:KTLX":'--',  "mult=1.15,noise=0.75:KTLX":':', 
    "adapt=0.80:KTLX":'-', "adapt=0.80,noise=0.5:KTLX":'--', "adapt=0.80,noise=0.75:KTLX":':', 
    "relax=0.50:KTLX":'-', "relax=0.50,noise=0.5:KTLX":'--', "relax=0.50,noise=0.75:KTLX":':',
    "no-infl:KTLX":'-'
}

radar_id_1km_goshen = [ "KCYS" ] #, "KFTG", "05XP" ]
experiments_1km_goshen = [ "1kmf-sndr0h=25km", "1kmf-zs25-no-05XP", "1kmf-zs25-no-mm", "1kmf-zs25-no-mm-05XP", "1kmf-z-no-v2", "1kmf-z-no-snd" ] #, "1kmf-prtrgn=1", "1kmf-newbc", "1kmf-snd-no-w", "1kmf-no-05XP", "1kmf-no-mm-05XP" ]
colors_1km_goshen = {
    "1kmf-sndr0h=25km:KCYS":'k',
    "1kmf-sndr0h=25km:KFTG":'k',
    "1kmf-sndr0h=25km:05XP":'k',
    "1kmf-zs25-no-05XP:KCYS":'r',
    "1kmf-zs25-no-05XP:KFTG":'r',
    "1kmf-zs25-no-05XP:05XP":'r',
    "1kmf-zs25-no-mm:KCYS":'g',
    "1kmf-zs25-no-mm:KFTG":'g',
    "1kmf-zs25-no-mm:05XP":'g',
    "1kmf-zs25-no-mm-05XP:KCYS":'b',
    "1kmf-zs25-no-mm-05XP:KFTG":'b',
    "1kmf-zs25-no-mm-05XP:05XP":'b',
    "1kmf-z-no-v2:KCYS":'m',
    "1kmf-z-no-v2:KFTG":'m',
    "1kmf-z-no-v2:05XP":'m',
    "1kmf-z-no-snd:KCYS":'c',
    "1kmf-z-no-snd:KFTG":'c',
    "1kmf-z-no-snd:05XP":'c',
}
styles_1km_goshen = {
    "1km-control-mod-05XP:KCYS":'-',
    "1km-control-mod-05XP:KFTG":'--',
    "1km-control-mod-05XP:05XP":':',
    "1km-control-no-mm:KCYS":'-',
    "1km-control-no-mm:KFTG":'--',
    "1km-control-no-mm:05XP":':',
    "1km-control-mm:KCYS":'-',
    "1km-control-mm:KFTG":'--',
    "1km-control-mm:05XP":':',
}
exp_names_1km_goshen = { 
    "1kmf-control":"Control Boundary Conditions",
    "1kmf-prtrgn=1":"Initial Perturbations in Storm Only", 
    "1kmf-newbc":"Control", 
    "1kmf-r0h=4km":"$r_{0h}$ = 4 km",
    "1kmf-snd-no-w":"Sounding obs do not update $w$",
    "1kmf-no-05XP":"No MWR-05XP data",
    "1kmf-no-mm-05XP":"No MM or MWR-05XP data",
    "1kmf-bc7dBZ,5ms":r"BC: $\sigma_Z$ = 7 dBZ, $\sigma_{v_r}$ = 5 m s$^{-1}$",
    "1kmf-bcmult=1.03":r"BC: Multiplicative $\alpha$ = 1.03",
    "1kmf-sndr0h=25km":r"Control (All VORTEX2 Obs)",
    "1kmf-zs25-no-05XP":"No MWR-05XP",
    "1kmf-zs25-no-mm":"No MM",
    "1kmf-zs25-no-mm-05XP":"No MWR-05XP or MM",
    "1kmf-z-no-v2":"No VORTEX2 Obs",
    "1kmf-z-no-snd":"No Soundings"
}

radar_id_3km_goshen = [ "KCYS", "KFTG", "KRIW" ]
experiments_3km_goshen = [ "3km-fixed-radar", "3kmf-7dBZ,5ms", "3kmf-r0h=12km", "3kmf-mult=1.03" ] #"3kmf-r0h=18km", "3kmf-pr0h=16km", "3kmf-posnegpt", "3kmf-pr0h=16km,r0h=12km", "3kmf-n0r=2e6" ]
colors_3km_goshen = {
    "3kmf-mult=1.03:KCYS":'r',
    "3kmf-mult=1.03:KFTG":'r',
    "3kmf-mult=1.03:KRIW":'r',
    "3kmf-r0h=12km:KCYS":'g',
    "3kmf-r0h=12km:KFTG":'g',
    "3kmf-r0h=12km:KRIW":'g',
    "3kmf-r0h=18km:KCYS":'r',
    "3kmf-r0h=18km:KFTG":'r',
    "3kmf-r0h=18km:KRIW":'r',
    "3kmf-pr0h=16km,r0h=12km:KCYS":'b',
    "3kmf-pr0h=16km,r0h=12km:KFTG":'b',
    "3kmf-pr0h=16km,r0h=12km:KRIW":'b',
    "3kmf-7dBZ,5ms:KCYS":'m',
    "3kmf-7dBZ,5ms:KFTG":'m',
    "3kmf-7dBZ,5ms:KRIW":'m',
    "3kmf-n0r=2e6:KCYS":'c',
    "3kmf-n0r=2e6:KFTG":'c',
    "3kmf-n0r=2e6:KRIW":'c',
    "3kmf-pr0h=16km:KCYS":'#999999',
    "3kmf-pr0h=16km:KFTG":'#999999',
    "3kmf-pr0h=16km:KRIW":'#999999',
    "3kmf-posnegpt:KCYS":'#ff6600',
    "3kmf-posnegpt:KFTG":'#ff6600',
    "3kmf-posnegpt:KRIW":'#ff6600',
    "3km-fixed-radar:KCYS":'k',
    "3km-fixed-radar:KFTG":'k',
    "3km-fixed-radar:KRIW":'k',
}
styles_3km_goshen = {
    "3km-control:KCYS":'-',
    "3km-control:KFTG":'--',
    "3km-control:KRIW":':',
    "3km-control-adapt=0.80:KCYS":'-',
    "3km-control-adapt=0.80:KFTG":'--',
    "3km-control-adapt=0.80:KRIW":':',
    "3km-control-adapt=1.00:KCYS":'-',
    "3km-control-adapt=1.00:KFTG":'--',
    "3km-control-adapt=1.00:KRIW":':',
}
exp_names_3km_goshen = { 
    "3kmf-control":r"Control (New)", 
    "3km-control-adapt=0.80":r"RTPS $\alpha$ = 0.80", 
    "3km-control-adapt=1.00":r"RTPS $\alpha$ = 1.00", 
    "3km-control-r0h=12km":r"$r_{0h}$ = 12 km",
    "3kmf-r0h=12km":r"$r_{0h}$ = 12 km",
    "3kmf-r0h=18km":r"$r_{0h}$ = 18 km",
    "3kmf-pr0h=16km":r"$r_{0h,i}$ = 16 km",
    "3kmf-pr0h=16km,r0h=12km":r"$r_{0h,i}$ = 16 km, $r_{0h}$ = 12 km",
    "3kmf-7dBZ,5ms":r"$\sigma_Z$ = 7 dBZ, $\sigma_{v_r}$ = 5 m s$^{-1}$",
    "3km-n0r=8e5":r"$N_{0r}$ = 8 $\times$ 10$^5$ m$^{-4}$",
    "3kmf-n0r=2e6":r"$N_{0r}$ = 2 $\times$ 10$^6$ m$^{-4}$",
    "3km-alladapt":"RTPS only",
    "3km-fixed-radar":"Control",
    "3kmf-posnegpt":r"+ and - $\theta$ perturbations",
    "3kmf-mult=1.03":r"Multiplicative $\alpha$ = 1.03 domain-wide",
}

def _loadRMSDFile(file_name):
    try:
        return [ float(v) for v in open(file_name, 'r').read().strip().split() ]
    except IOError:
        return tuple([np.nan, np.nan])

def _loadSpreadFile(file_name):
    try:
        return [ float(v) for v in open(file_name, 'r').read().strip().split() ]
    except IOError:
        return tuple([np.nan, np.nan])

def _loadconsistencyparameterfile(file_name):
    try:
        return [ float(v) for v in open(file_name, 'r').readline().strip().split() ]
    except IOError:
        return tuple([np.nan, np.nan])

def plotSubplots(times, data, legend_loc, y_label, y_lim, colors, styles, exp_names, title, file_name):
    pylab.figure(figsize=(10, 8))
    pylab.subplots_adjust(left=0.075, right=0.95, top=0.925, bottom=0.1, wspace=0.175, hspace=0.275)

    radars = sorted(list(set([ name.split(":")[-1] for name in data.keys() ])))

    n_rows = 1
    n_cols = (len(radars) + 1) / 2

    lines = {}

    for idx, radar in enumerate(radars):
        all_good = []

        pylab.subplot(n_rows, n_cols, idx + 1)
        for exp in sorted([ e for e in data.keys() if e.split(":")[-1] == radar ]):
            good_idxs = np.where(~np.isnan(data[exp]))[0]
            if len(good_idxs) > 0:
                name = exp.split(":")[0]
                line = pylab.plot(times[good_idxs], data[exp][good_idxs], color=colors[exp], label=exp_names[name])
                lines[exp_names[name]] = line

            all_good.append(good_idxs)
        all_good_idxs = np.unique1d(np.concatenate(tuple(all_good)))
        pylab.plot([times.min(), times.max()], [0, 0], color='k', linestyle=':')
        pylab.axvline(14400, color='k', linestyle=':')

        pylab.xlabel(r"Time (UTC)", size='xx-large')
        pylab.ylabel(y_label, size='xx-large')
        pylab.xlim((times.min(), times.max()))
        pylab.ylim(y_lim)
        unique_times = np.sort(np.unique1d(times))
        pylab.xticks(unique_times[::2], [ (datetime(2009, 6, 5, 18, 0, 0) + timedelta(seconds=int(t))).strftime("%H%M") for t in unique_times ][::2], rotation=30, size='xx-large')
        pylab.yticks(size='xx-large')

#       pylab.title(radar)

    labels, line_objs = zip(*lines.items())
#   pylab.gcf().legend(line_objs, labels, 'lower right', prop={'size':'medium'})
    pylab.legend(line_objs, labels, loc=legend_loc, prop={'size':'medium'})

    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def plot(times, rms_difference, legend_loc, y_label, y_lim, colors, styles, title, file_name):
#   exp_names = { "1km-control-mod-05XP":"MM + MWR05XP", "1km-control-no-mm":"No MM", "1km-control-mm":"MM" }
    exp_names = { "3km-control":r"5 dBZ, 3 m s$^{-1}$", "3km-control-adapt=0.80":r"RTPS $\alpha$ = 0.80", "3km-control-adapt=1.00":r"RTPS $\alpha$ = 1.00", 
        "3km-control-r0h=12km":r"$r_{0h}$ = 12 km", "3km-control-7dBZ,5ms":r'$\sigma_Z$ = 7 dBZ, $\sigma_{v_r}$ = 5 m s$^{-1}$' }
    pylab.figure()
    pylab.axes((0.1, 0.125, 0.85, 0.8))

    all_good = [] 

    for exp_name in sorted(rms_difference.keys()):
        good_idxs = np.where(~np.isnan(rms_difference[exp_name]))[0]

        name, radar= exp_name.split(':')
        if len(good_idxs) > 0:
            pylab.plot(times[good_idxs], rms_difference[exp_name][good_idxs], color=colors[exp_name], linestyle=styles[exp_name], label="%s (%s)" % (exp_names[name], radar))

        all_good.append(good_idxs)

    all_good_idxs = np.unique1d(np.concatenate(tuple(all_good)))

    pylab.plot([times.min(), times.max()], [0, 0], color='k', linestyle=':')

    pylab.xlabel(r"Time (UTC)", size='large')
    pylab.ylabel(y_label, size='large')
    pylab.xlim((times.min(), times.max()))
    pylab.ylim(y_lim)
    pylab.xticks(times[all_good_idxs], [ (datetime(2009, 6, 5, 18, 0, 0) + timedelta(seconds=int(t))).strftime("%H%M") for t in times[all_good_idxs] ], size='large', rotation=30)
    pylab.yticks(size='large')
    pylab.legend(loc=legend_loc, prop={'size':'small'})
    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def main():
    base_path = "/caps2/tsupinie/"
#   base_path = "/data6/tsupinie/goshen/"

    tag = "1km"

    if tag == "cov_infl":
        radar_id = radar_id_cov_infl
        experiments = experiments_cov_infl
        colors = colors_cov_infl
        styles = styles_cov_infl

        t_ens_start = 0
        t_ens_end = 3600
        t_ens_step = 300
    elif tag == "1km":
        radar_id = radar_id_1km_goshen
        experiments = experiments_1km_goshen
        colors = colors_1km_goshen
        styles = styles_1km_goshen
        exp_names = exp_names_1km_goshen

        t_ens_start = 10800
        t_ens_end = 18000
        t_ens_step = 300
    elif tag == "3km":
        radar_id = radar_id_3km_goshen
        experiments = experiments_3km_goshen
        colors = colors_3km_goshen
        styles = styles_3km_goshen
        exp_names = exp_names_3km_goshen

        t_ens_start = 10800
        t_ens_end = 18000
        t_ens_step = 300

    rmsd = {'vr':{}, 'ref':{}}
    spread = {'vr':{}, 'ref':{}}
    consistency = {'vr':{}, 'ref':{}}

    rmsd_sawtooth = {'vr':{}, 'ref':{}}
    spread_sawtooth = {'vr':{}, 'ref':{}}

    for rid in radar_id:
        for exp_name in experiments:
            key = "%s:%s" % (exp_name, rid)

            for quant in ['vr', 'ref']:
                rmsd[quant][key] = []
                spread[quant][key] = []
                consistency[quant][key] = []

                rmsd_sawtooth[quant][key] = []
                spread_sawtooth[quant][key] = []

            for t_ens in xrange(t_ens_start + t_ens_step, t_ens_end + t_ens_step, t_ens_step):
#               if exp_name[:4] == "mult":
#                   rmsd_file_name = "%s/%s/%srmsdbg%06d" % (base_path, exp_name, rid, t_ens)
#               else:
#                   rmsd_file_name = "%s/%s/%srmsdbg%06d" % (base_path, exp_name, rid, t_ens)

                rmsd_file_name = "%s/%s/%srmsdbg%06d" % (base_path, exp_name, rid, t_ens)
                vr_rmsd, ref_rmsd = _loadRMSDFile(rmsd_file_name)
                rmsd['vr'][key].append(vr_rmsd)
                rmsd['ref'][key].append(ref_rmsd)

                spread_file_name = "%s/%s/%sspreadfcs%06d" % (base_path, exp_name, rid, t_ens)
                vr_spread, ref_spread = _loadSpreadFile(spread_file_name)
                spread['vr'][key].append(vr_spread)
                spread['ref'][key].append(ref_spread)

                consistency_file_name = "%s/%s/%sinnov%06d" % (base_path, exp_name, rid, t_ens)
                vr_consistency, ref_consistency = _loadConsistencyParameterFile(consistency_file_name)
                consistency['vr'][key].append(vr_consistency)
                consistency['ref'][key].append(ref_consistency)

            for t_ens in xrange(t_ens_start, t_ens_end + t_ens_step, t_ens_step):
                rmsd_file_name = "%s/%s/%srmsdbg%06d" % (base_path, exp_name, rid, t_ens)
                vr_rmsd, ref_rmsd = _loadRMSDFile(rmsd_file_name)
                rmsd_sawtooth['vr'][key].append(vr_rmsd)
                rmsd_sawtooth['ref'][key].append(ref_rmsd)
                
                rmsd_file_name = "%s/%s/%srmsdan%06d" % (base_path, exp_name, rid, t_ens)
                vr_rmsd, ref_rmsd = _loadRMSDFile(rmsd_file_name)
                rmsd_sawtooth['vr'][key].append(vr_rmsd)
                rmsd_sawtooth['ref'][key].append(ref_rmsd)

                spread_file_name = "%s/%s/%sspreadfcs%06d" % (base_path, exp_name, rid, t_ens)
                vr_spread, ref_spread = _loadRMSDFile(spread_file_name)
                spread_sawtooth['vr'][key].append(vr_spread)
                spread_sawtooth['ref'][key].append(ref_spread)
                
                spread_file_name = "%s/%s/%sspreadana%06d" % (base_path, exp_name, rid, t_ens)
                vr_spread, ref_spread = _loadRMSDFile(spread_file_name)
                spread_sawtooth['vr'][key].append(vr_spread)
                spread_sawtooth['ref'][key].append(ref_spread)

            for quant in ['vr', 'ref']:
                rmsd[quant][key] = np.array(rmsd[quant][key])
                spread[quant][key] = np.array(spread[quant][key])
                consistency[quant][key] = np.log10(np.array(consistency[quant][key]))

                rmsd_sawtooth[quant][key] = np.array(rmsd_sawtooth[quant][key])
                spread_sawtooth[quant][key] = np.array(spread_sawtooth[quant][key])

                if quant == 'ref' and rid == "05XP":
                    rmsd[quant][key] = np.nan * np.zeros(rmsd[quant][key].shape)
                    spread[quant][key] = np.nan * np.zeros(spread[quant][key].shape)
                    consistency[quant][key] = np.nan * np.zeros(consistency[quant][key].shape)

                    rmsd_sawtooth[quant][key] = np.nan * np.zeros(rmsd_sawtooth[quant][key].shape)
                    spread_sawtooth[quant][key] = np.nan * np.zeros(spread_sawtooth[quant][key].shape)

    times = np.arange(t_ens_start + t_ens_step, t_ens_end + t_ens_step, t_ens_step)

#   plot(times, rmsd['ref'], 1, "RMS Difference (dBZ)", (0, 25), colors, styles, "RMS Analysis Innovation for Reflectivity", "rmsd_ref_%s.png" % tag)
#   plot(times, spread['ref'], 2, "Spread (dBZ)", (0, 25), colors, styles, "Ensemble Analysis Spread in Reflectivity", "spread_ref_%s.png" % tag)
    plotSubplots(times, consistency['ref'], 1, "Log Consistency Parameter (unitless)", (-1, 1), colors, styles, exp_names, "Log Consistency Parameter for Reflectivity", "consistency_ref_%s.png" % tag)

#   plot(times, rmsd['vr'], 1, "RMS Difference (m s$^{-1}$)", (0, 10), colors, styles, "RMS Analysis Innovation for Radial Velocity", "rmsd_vr_%s.png" % tag)
#   plot(times, spread['vr'], 2, "Spread (m s$^{-1}$)", (0, 10), colors, styles, "Ensemble Analysis Spread in Radial Velocity", "spread_vr_%s.png" % tag)
    plotSubplots(times, consistency['vr'], 3, "Log Consistency Parameter (unitless)", (-1, 1), colors, styles, exp_names, "Log Consistency Parameter for Radial Velocity", "consistency_vr_%s.png" % tag)

    times = np.arange(t_ens_start, t_ens_end + t_ens_step, t_ens_step)
    plotSubplots(times[:, np.newaxis].repeat(2, axis=1).flatten('C'), rmsd_sawtooth['ref'], 4, "RMS Difference (dBZ)", (0, 20), colors, styles, exp_names, "RMS Innovation for Reflectivity", "rmsd_sawtooth_ref_%s.png" % tag)
    plotSubplots(times[:, np.newaxis].repeat(2, axis=1).flatten('C'), rmsd_sawtooth['vr'], 2, "RMS Difference (m s$^{-1}$)", (0, 10), colors, styles, exp_names, "RMS Innovation for Radial Velocity", "rmsd_sawtooth_vr_%s.png" % tag)

    plotSubplots(times[:, np.newaxis].repeat(2, axis=1).flatten('C'), spread_sawtooth['ref'], 1, "Spread (dBZ)", (0, 20), colors, styles, exp_names, "Spread for Reflectivity", "spread_sawtooth_ref_%s.png" % tag)
    plotSubplots(times[:, np.newaxis].repeat(2, axis=1).flatten('C'), spread_sawtooth['vr'], 1, "Spread (m s$^{-1}$)", (0, 10), colors, styles, exp_names, "Spread for Radial Velocity", "spread_sawtooth_vr_%s.png" % tag)
    return

if __name__ == "__main__":
    main()
