
import cPickle
from datetime import datetime, timedelta
import warnings

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap

from legacy import setupMapProjection, goshen_1km_proj, goshen_1km_gs, inflow_stations
from util import loadObs

def partitionObs(obs, base_epoch):
    times = sorted(inflow_stations.keys())
    regions = inflow_stations[times[0]].keys()
    partitions = dict([ (r, np.empty((0,), dtype=obs.dtype)) for r in regions ])

    for t_ens in times:
        t_obs = obs[np.where(obs['time'] == base_epoch + t_ens)]

        ids = np.unique1d(t_obs['id'])
        for id in ids:
            for region in regions:
                if id in inflow_stations[t_ens][region]:
                    ob = t_obs[np.where(t_obs['id'] == id)]
                    partitions[region] = np.append(partitions[region], ob)

    return partitions

def compile(data, var_name, exp_name=None, region=None):
    compiled_data = []

    for name, exps in data.iteritems():
        for var, variables in exps.iteritems():
            for reg, regions in variables.iteritems():
                if (exp_name is None or name == exp_name) and (var == var_name) and (region is None or reg == region):
                    compiled_data.extend(regions)

    return np.array(compiled_data)

def reliability(alpha_bar, beta_bar, high_outlier, low_outlier):
    gs = alpha_bar + beta_bar
    os = beta_bar / (alpha_bar + beta_bar)

    os[0] = low_outlier
    os[-1] = high_outlier

    gs[0] = beta_bar[0] / os[0]
    gs[-1] = alpha_bar[-1] / (1 - os[-1])
    gs = np.where(np.isnan(gs), 0, gs)

    prob = np.arange(alpha_bar.shape[0], dtype=float) / alpha_bar.shape[0]

    return (gs * (os - prob) ** 2).sum()

def uncertainty(obs):
    obs_sort = np.sort(obs)
    probs = np.arange(1, obs.shape[0], dtype=float) / obs.shape[0]

    return (probs * (1 - probs) * (obs_sort[1:] - obs_sort[:-1])).sum()

def windDirSpd2UV(wind_dir, wind_spd):
    u = -wind_spd * np.sin(np.radians(wind_dir))
    v = -wind_spd * np.cos(np.radians(wind_dir))
    return u, v

def computeUncertainty(var_order, region_order):
    base_time = datetime(2009, 6, 5, 18, 0, 0)
    epoch = datetime(1970, 1, 1, 0, 0, 0)
    base_epoch = (base_time - epoch).total_seconds()
    times = np.arange(14700, 18300, 300)

    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs)
    map = Basemap(**proj)

    obs_file_names = ['psu_straka_mesonet.pkl', 'ttu_sticknet.pkl', 'asos.pkl', 'soundings_clip.pkl']
    obs = loadObs(obs_file_names, [ base_time + timedelta(seconds=int(t)) for t in times ], map, (goshen_1km_proj['width'], goshen_1km_proj['height']), sounding_obs=['soundings_clip.pkl'])


    obs_part = partitionObs(obs, base_epoch)
    u, v = windDirSpd2UV(obs['wind_dir'], obs['wind_spd'])
    u_part, v_part = {}, {}
    for region, reg_obs in obs_part.iteritems():
        u_part[region], v_part[region] = windDirSpd2UV(obs_part[region]['wind_dir'], obs_part[region]['wind_spd'])

    all_uncert = []

    def row(all_obs, part_obs, name, units):
        uncert = uncertainty(all_obs)
        row_uncert = [ uncert ]

        row = "\t" + r"%s (%s) & %.2f" % (name, units, uncert)
        for region in region_order:
            uncert = uncertainty(part_obs[region])
            row += " & %.2f" % uncert
            row_uncert.append(uncert)
        print row + r"\\"
        print "\t" + r"\hline"
        return np.array(row_uncert)

    all_uncert.append(row(obs['temp'], dict([ (key, val['temp'] ) for key, val in obs_part.iteritems()]), '$T$', r'$^{\circ}$F'))
    all_uncert.append(row(obs['dewp'], dict([ (key, val['dewp'] ) for key, val in obs_part.iteritems()]), '$T_d$', r'$^{\circ}$F'))
    all_uncert.append(row(u, u_part, '$u$', r'm s$^{-1}$'))
    all_uncert.append(row(v, v_part, '$v$', r'm s$^{-1}$'))
    return np.array(all_uncert)

def computeReliability(var_order, exp_order, region_order):
    all_alphas, all_betas, all_highs, all_lows = {}, {}, {}, {}

    for exp_name in exp_order:
        all_alphas[exp_name], all_betas[exp_name], all_highs[exp_name], all_lows[exp_name] = cPickle.load(open("%s_crps_pieces.pkl" % exp_name, 'r'))

    units = { 't':r"$^{\circ}$F", 'td':r"$^{\circ}$F", 'u':r"m s$^{-1}$", 'v':r"m s$^{-1}$" }
    var_names = { 't':'T', 'td':"T_d", 'u':'u', 'v':'v' }
    exp_labels = { 'no-mm':"No MM", 'mm':"MM", 'mod-05XP':"MM+MWR05XP" }

    all_reli = []
    for var in var_order:
        print "\t" + r"\multirow{3}{*}{$%s$ (%s)}" % (var_names[var], units[var])
        var_reli = []
        for exp_name in exp_order:
            compiled_alphas, compiled_betas, compiled_highs, compiled_lows = [ compile(p, var, exp_name=exp_name).mean(axis=0) for p in [ all_alphas, all_betas, all_highs, all_lows ] ]
            reli = reliability(compiled_alphas, compiled_betas, compiled_highs, compiled_lows)
            row_reli = [ reli ]
            row = "\t& %s & %.2f " % (exp_labels[exp_name], reli)
            for region in region_order:
                compiled_alphas, compiled_betas, compiled_highs, compiled_lows = [ compile(p, var, exp_name=exp_name, region=region).mean(axis=0) for p in [ all_alphas, all_betas, all_highs, all_lows ] ]
                reli = reliability(compiled_alphas, compiled_betas, compiled_highs, compiled_lows)
                row += "& %.2f " % reli

                row_reli.append(reli)
            row += r"\\"
            print row
            var_reli.append(np.array(row_reli))

        print "\t" + r"\hline"
        all_reli.append(np.array(var_reli))
    return np.array(all_reli)

def computeCRPS(var_order, exp_order, region_order):
    crps_data = {}

    for name in exp_order:
        crps_data[name] = cPickle.load(open("%s_crps_breakdown.pkl" % name, 'r'))

    units = { 't':r"$^{\circ}$F", 'td':r"$^{\circ}$F", 'u':r"m s$^{-1}$", 'v':r"m s$^{-1}$" }
    var_names = { 't':'T', 'td':"T_d", 'u':'u', 'v':'v' }
    exp_labels = { 'no-mm':"No MM", 'mm':"MM", 'mod-05XP':"MM+MWR05XP" }

    all_crps = []
    for var in var_order:
        print "\t" + r"\multirow{3}{*}{$%s$ (%s)}" % (var_names[var], units[var])

        var_crps = []
        for name in exp_order:
            crps = compile(crps_data, var, exp_name=name).mean(axis=0)
            row_crps = [ crps]
            
            row = "\t& %s & %.2f " % (exp_labels[name], crps)
            for region in region_order:
                crps = compile(crps_data, var, exp_name=name, region=region).mean(axis=0)
                row += "& %.2f " % crps
                row_crps.append(crps)
            row += r"\\"
            print row
            var_crps.append(np.array(row_crps))

        print "\t" + r"\hline"
        all_crps.append(np.array(var_crps))
    return np.array(all_crps)

def main():
    units = { 't':r"$^{\circ}$F", 'td':r"$^{\circ}$F", 'u':r"m s$^{-1}$", 'v':r"m s$^{-1}$" }
    var_names = { 't':'T', 'td':"T_d", 'u':'u', 'v':'v' }
    exp_labels = { 'no-mm':"No MM", 'mm':"MM", 'mod-05XP':"MM+MWR05XP" }

    exp_order = [ 'no-mm', 'mm', 'mod-05XP' ]
    region_order = ['sounding', 'inflow', 'outflow']
    var_order = ['t', 'td', 'u', 'v']
    print " *** LaTex Table for uncertainty *** "
    uncert = computeUncertainty(var_order, region_order)
    print " *** LaTex Table for reliability *** "
    reli = computeReliability(var_order, exp_order, region_order)
    print " *** LaTex Table for CRPS *** "
    crps = computeCRPS(var_order, exp_order, region_order)

    reli = np.transpose(reli, (0, 2, 1))
    crps = np.transpose(crps, (0, 2, 1))

    print uncert.shape
    print reli.shape
    print crps.shape

    resol = reli + uncert[:, :, np.newaxis] - crps
    resol = np.transpose(resol, (0, 2, 1))

    print " *** LaTex Table for resolution *** " 
    for idx, var in enumerate(var_order):
        print "\t" + r"\multirow{3}{*}{$%s$ (%s)}" % (var_names[var], units[var])
        for jdy, exp_name in enumerate(exp_order):
            row = "\t& %s & %.2f " % (exp_labels[exp_name], resol[idx, jdy, 0])
            for kdz, region in enumerate(region_order):
                row += "& %.2f " % resol[idx, jdy, kdz + 1]

            row += r"\\"
            print row

        print "\t" + r"\hline"
    return

if __name__ == "__main__":
    main()
