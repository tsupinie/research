
from util import setupMapProjection, decompressVariable, loadObs, goshen_1km_proj, goshen_1km_gs, inflow_stations
from dataload import loadEnsemble

import numpy as np

import Nio as nio

import matplotlib
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap

import cPickle, glob, os
from datetime import datetime, timedelta

def heaviside(x):
    if x >= 0:
        return 1
    else:
        return 0

def CRPS(ens_values, observation):
    sort_ens_values = np.sort(ens_values)

    bin = binRank(ens_values, observation)
    prob = np.arange(ens_values.shape[0] + 1, dtype=float) / len(ens_values)
    alphas = np.zeros((ens_values.shape[0] + 1,), dtype=float)
    betas = np.zeros((ens_values.shape[0] + 1,), dtype=float)

    if bin != 0:
        alphas[1:bin] = sort_ens_values[1:bin] - sort_ens_values[:(bin - 1)]
        alphas[bin] = observation - sort_ens_values[bin - 1]

    if bin != len(ens_values):
        betas[(bin + 1):-1] = sort_ens_values[(bin + 1):] - sort_ens_values[bin:-1]
        betas[bin] = sort_ens_values[bin] - observation

    return (alphas * prob ** 2 + betas * (1 - prob) ** 2).sum(), alphas, betas

def binRank(ens_values, observation):
    sort_ens_values = np.sort(ens_values)
    closest_idx = np.argmin(np.abs(sort_ens_values - observation))
    if observation > sort_ens_values[closest_idx]:
        closest_idx += 1
    return closest_idx

def tempFromPt(**kwargs):
    temp = kwargs['pt'] * (kwargs['p'] / 100000.) ** (2. / 7.)
    return (temp - 273.15) * 9. / 5. + 32

def dewpFromQv(**kwargs):
    e = (kwargs['p'] * kwargs['qv']) / (kwargs['qv'] + 0.622)
    dewp = 1. / (1 / 273.15 - 461.5 / 2.5e6 * np.log(e / 611.))
    return (dewp - 273.15) * 9. / 5. + 32

def uFromU(**kwargs):
    return kwargs['u']

def vFromV(**kwargs):
    return kwargs['v']

def uFromWind(**kwargs):
    return -kwargs['wind_spd'] * np.sin(kwargs['wind_dir'] * np.pi / 180.)

def vFromWind(**kwargs):
    return -kwargs['wind_spd'] * np.cos(kwargs['wind_dir'] * np.pi / 180.)

def tempFromT(**kwargs):
    return kwargs['temp']

def dewpFromTd(**kwargs):
    return kwargs['dewp']

def plotCDFs(ens_members, observation, title, file_name):
    margin_factor = 1.1

    x_lb = min(ens_members.min(), observation)
    x_ub = max(ens_members.max(), observation)

    x_ctr = (x_lb + x_ub) / 2
    
    x_lb = x_ctr - margin_factor * (x_ctr - x_lb)
    x_ub = x_ctr + margin_factor * (x_ub - x_ctr)

    pylab.figure(figsize=(10, 6))
    ens_members_expanded = np.vstack((np.insert(ens_members, 0, x_lb), np.append(ens_members, x_ub))).flatten('F')

    fcst_cdf = np.arange(ens_members.shape[0] + 1, dtype=float) / ens_members.shape[0]
    fcst_cdf_expanded = np.vstack((fcst_cdf, fcst_cdf)).flatten('F')

    pylab.plot(ens_members_expanded, fcst_cdf_expanded, 'b-')
    pylab.plot([x_lb, observation, observation, x_ub], [0, 0, 1, 1], 'k-')

    pylab.xlim(x_lb, x_ub)
    pylab.ylim(-0.1, 1.1)

    pylab.title(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def createVerificationGraphs(alphas, betas, high_outliers, low_outliers, rank_histogram, n_obs, suffix, exp_name):
    if suffix != "":
        suffix = "_%s" % suffix

    gs = alphas + betas
    os = betas / (alphas + betas)

    os[0] = low_outliers
    os[-1] = high_outliers

    gs[0] = betas[0] / os[0]
    gs[-1] = alphas[-1] / (1 - os[-1])

    gs = np.where(np.isnan(gs), np.zeros(gs.shape), gs)

    prob = np.arange(alphas.shape[0], dtype=np.float) / (alphas.shape[0] - 1)
    prob_expanded = np.vstack((prob, prob)).flatten('F')[1:]
    reliability = (gs * (os - prob) ** 2).sum()

    cum_gs = np.add.accumulate(gs)
    cum_gs_expanded = np.vstack((cum_gs, cum_gs)).flatten('F')[:-1]
    pylab.clf()
    pylab.plot(cum_gs_expanded, prob_expanded)
    pylab.savefig("images-%s/crps_cum_bin_width%s.png" % (exp_name, suffix))

    pylab.clf()
    pylab.plot(prob, os)
    pylab.plot([0.0, 1.0], [0.0, 1.0], 'k', lw=0.5)
    pylab.text(0.975, 0.025, "Reliability: %.3f" % reliability, style='italic', size='large', weight='bold', transform=pylab.gca().transAxes, ha='right', bbox={'facecolor':'#ffffff'})
    pylab.xlabel("Normalized Rank")
    pylab.ylabel("Observed Frequency")
    pylab.savefig("images-%s/crps_reliability%s.png" % (exp_name, suffix))

    pylab.clf()
    pylab.bar(np.arange(rank_histogram.shape[0]), rank_histogram)
    pylab.xlim(0, rank_histogram.shape[0])
#   pylab.ylim(0.0, 1.0)
#   pylab.text(0.5, 0.8, "Number of\nobservations: %d" % n_obs, style='italic', size='large', weight='bold', transform=pylab.gca().transAxes, ha='center', va='center', bbox={'facecolor':'#ffffff'})
    pylab.xlabel("Rank", size='large')
    pylab.ylabel("Normalized Bin Count (n = %d)" % n_obs, size='large')
    pylab.xticks(size='large')
    pylab.yticks(size='large')
    pylab.savefig("images-%s/crps_rankhist%s.png" % (exp_name, suffix))
    pylab.close()
    return

def rearrangeSoundingObs(sounding_obs):
    dtype = [('id', np.dtype('|S4')), ('obtype', np.dtype('|S8')), ('time', np.dtype('float64')), ('latitude', np.dtype('float64')), ('longitude', np.dtype('float64')), ('elevation', np.dtype('float64')), 
                ('temp', np.dtype('float64')), ('dewp', np.dtype('float64')), ('pres', np.dtype('float64')), ('wind_dir', np.dtype('float64')), ('wind_spd', np.dtype('float64'))]

    rearr_sounding_obs = []
    for sounding in sounding_obs:
        release_time = datetime.strptime(sounding['release_time'], "%Y, %m, %d, %H:%M:%S")
        release_epoch = (release_time - datetime(1970, 1, 1, 0, 0, 0)).total_seconds()
        for time, lat, lon, hght, temp, dewp, pres, wdir, wspd in zip(*[ sounding[k] for k in ['time', 'latitude', 'longitude', 'altitude', 'temperature', 'dewpoint', 'pressure', 'wind_dir', 'wind_spd']]):
            if temp != 999.0 and dewp != 999.0 and pres != 9999.0 and wdir != 999.0 and wspd != 999.0:
               rearr_sounding_obs.append((sounding['release_site'][:4], " " * 8, release_epoch + time, lat, lon, hght, temp, dewp, pres, wdir, wspd))

    return np.array(rearr_sounding_obs, dtype=dtype)

def main():
    base_time = datetime(2009, 6, 5, 18, 0, 0)
    epoch = datetime(1970, 1, 1, 0, 0, 0)
    times_seconds = range(14700, 18300, 300)
    times = [ base_time + timedelta(seconds=t) for t in times_seconds ]

    n_ensemble_members = 40
    exp_name = "zupdtpt"

    #
    # Set up the basemap grid
    #
    proj = setupMapProjection(goshen_1km_proj, goshen_1km_gs)
    map = Basemap(**proj)

    #
    # Load and thin all the observed data
    #
    obs_file_names = ['psu_straka_mesonet.pkl', 'ttu_sticknet.pkl', 'asos.pkl', 'soundings_clip.pkl']
    all_obs = loadObs(obs_file_names, times, map, (goshen_1km_proj['width'], goshen_1km_proj['height']), sounding_obs=['soundings_clip.pkl'])
    print all_obs.shape[0]

    ob_first_char = np.array([ id[0] for id in list(all_obs['id']) ])

    num_psu_obs = len(np.where(ob_first_char == "P")[0])
    num_ttu_obs = len(np.where((ob_first_char == "1") | (ob_first_char == "2"))[0])
    num_asos_obs = len(np.where((ob_first_char == "K"))[0])
    num_sndg_obs = len(np.where(all_obs['obtype'] == "SNDG")[0])
    print "Number of NSSL MM obs used:", num_psu_obs
    print "Number of TTU Sticknet obs used:", num_ttu_obs
    print "Number of ASOS obs used:", num_asos_obs
    print "Number of sounding obs used:", num_sndg_obs

    all_times = [ datetime(1970, 1, 1, 0, 0, 0) + timedelta(seconds=t) for t in all_obs['time'] ]

    #
    # Convert the latitude and longitude observations to x and y on the grid.
    #
    obs_x, obs_y = map(all_obs['longitude'], all_obs['latitude'])
    obs_z = all_obs['pres'] * 100

    def getObsData(**kwargs): 
        obs = np.empty(kwargs['pt'].shape, dtype=[('u', np.float32), ('v', np.float32), ('pt', np.float32), ('p', np.float32), ('qv', np.float32)])
        obs['u'] = kwargs['u']
        obs['v'] = kwargs['v']
        obs['pt'] = kwargs['pt']
        obs['p'] = kwargs['p']
        obs['qv'] = kwargs['qv']

        return obs

    obs_vars = ['u', 'v', 't', 'td']
    ens_funcs = { 'u':uFromU, 'v':vFromV, 't':tempFromPt, 'td':dewpFromQv }
    obs_funcs = { 'u':uFromWind, 'v':vFromWind, 't':tempFromT, 'td':dewpFromTd }

    avg_crps_values = { }
    all_crps_values = { }
    rank_histograms = { }
    all_alphas = { }
    all_betas = { }    
    high_outliers = { }
    low_outliers = { }

    for time_sec, time in zip(times_seconds, times):
#       files = glob.glob("/caps2/tsupinie/1kmf-%s/ena???.hdf%06d" % (exp_name, time_sec))

        time_idxs = np.where(all_obs['time'] == (time - epoch).total_seconds())

        #
        # Load all the ensemble members and interpolate them to the observation points.  Because of the design of my script, I'm
        # loading the all the members timestep-by-timestep, but there's no reason you can't load them all at once.  See the function
        # definition for the meaning of all the arguments.
        #
#       ens_obs, ens_members, ens_times = loadAndInterpolateEnsemble(files, ['u', 'v', 'pt', 'p', 'qv'], getObsData, "/caps2/tsupinie/1kmf-%s/ena001.hdfgrdbas" % exp_name, 
#           {'z':obs_z[time_idxs], 'y':obs_y[time_idxs], 'x':obs_x[time_idxs]}, agl=False, wrap=True, coords='pres')

        ens_obs = loadEnsemble("/caps2/tsupinie/1kmf-%s/" % exp_name, n_ensemble_members, [ time_sec ], (['u', 'v', 'pt', 'p', 'qv'], getObsData), { 'z':obs_z[time_idxs], 'y':obs_y[time_idxs], 'x':obs_x[time_idxs] }, agl=False, wrap=True, coords='pres')

#       print ens_obs

        #
        # All subsequent lines do the verification
        #
        for ob_var in obs_vars:
            time_crps_values = []

            ens_ob_var = ens_funcs[ob_var](**dict([ (n, ens_obs[n][:, 0]) for n in ens_obs.dtype.names ]))
            obs = obs_funcs[ob_var](**dict([ (n, all_obs[n][time_idxs]) for n in all_obs.dtype.names ]))

            if ob_var not in rank_histograms:
                rank_histograms[ob_var] = {}
                all_crps_values[ob_var] = {}
                all_alphas[ob_var] = {}
                all_betas[ob_var] = {}
                high_outliers[ob_var] = {}
                low_outliers[ob_var] = {}
               
                for region in [ 'inflow', 'outflow', 'sounding' ]:
                    rank_histograms[ob_var][region] = np.zeros((ens_obs.shape[0] + 1,), dtype=int)
                    all_crps_values[ob_var][region] = []
                    all_alphas[ob_var][region] = []
                    all_betas[ob_var][region] = []
                    high_outliers[ob_var][region] = []
                    low_outliers[ob_var][region] = []

            for idx in xrange(obs.shape[-1]):
                rank_idx = binRank(ens_ob_var[:, idx], obs[idx])
                crps, alphas, betas = CRPS(ens_ob_var[:, idx], obs[idx])
                high_outlier = heaviside(ens_ob_var[:, idx].max() - obs[idx])
                low_outlier = heaviside(ens_ob_var[:, idx].min() - obs[idx])

                for region in [ 'inflow', 'outflow', 'sounding' ]:
                    if region in inflow_stations[time_sec] and all_obs['id'][time_idxs][idx] in inflow_stations[time_sec][region]:
#                       plotCDFs(np.sort(ens_ob_var[:, idx]), obs[idx], "CDFs for Surface %s Observation %d and Forecast" % (ob_var, idx), "crps_cdf_sfc_%s_%03d.png" % (ob_var, idx))

                        rank_histograms[ob_var][region][rank_idx] += 1

                        all_crps_values[ob_var][region].append(crps)
                        all_alphas[ob_var][region].append(alphas)
                        all_betas[ob_var][region].append(betas)
                        high_outliers[ob_var][region].append(high_outlier)
                        low_outliers[ob_var][region].append(low_outlier)

                    elif region == "sounding" and all_obs['obtype'][time_idxs][idx] == "SNDG":
#                       plotCDFs(np.sort(ens_ob_var[:, idx]), obs[idx], "CDFs for Sounding %s Observation %d and Forecast" % (ob_var, idx), "crps_cdf_sndg_%s_%03d.png" % (ob_var, idx))

                        rank_histograms[ob_var][region][rank_idx] += 1

                        all_crps_values[ob_var][region].append(crps)
                        all_alphas[ob_var][region].append(alphas)
                        all_betas[ob_var][region].append(betas)
                        high_outliers[ob_var][region].append(high_outlier)
                        low_outliers[ob_var][region].append(low_outlier)

                time_crps_values.append(crps)

            try:
                avg_crps_values[ob_var].append(sum(time_crps_values) / len(time_crps_values))
            except KeyError:
                avg_crps_values[ob_var] = [ sum(time_crps_values) / len(time_crps_values) ]

    def dictmean(D):
        all_lists = []
        for val in D.itervalues(): all_lists.extend(val)
        return np.array(all_lists).mean(axis=0)

    def dictsum(D):
        all_lists = []
        for val in D.itervalues(): all_lists.append(val)
        return np.array(all_lists).sum(axis=0)

    def mean(L):
        return np.array(L).mean(axis=0)

    if not os.path.exists("images-%s" % exp_name):
        os.mkdir("images-%s" % exp_name, 0755)

    cPickle.dump(avg_crps_values, open("%s_crps.pkl" % exp_name, 'w'), -1)
    cPickle.dump(all_crps_values, open("%s_crps_breakdown.pkl" % exp_name, 'w'), -1)
    cPickle.dump((all_alphas, all_betas, high_outliers, low_outliers), open("%s_crps_pieces.pkl" % exp_name, 'w'), -1)

    for ob_var in obs_vars:
        total_obs = sum([ len(v) for v in high_outliers[ob_var].itervalues()  ])
        print total_obs
        createVerificationGraphs(dictmean(all_alphas[ob_var]), dictmean(all_betas[ob_var]), dictmean(high_outliers[ob_var]), dictmean(low_outliers[ob_var]), dictsum(rank_histograms[ob_var]).astype(float) / total_obs, total_obs, "%s" % ob_var, exp_name)

        for region in [ 'inflow', 'outflow',  'sounding' ]:
            suffix = "%s_%s" % (ob_var, region)
            region_obs = len(high_outliers[ob_var][region])
            createVerificationGraphs(mean(all_alphas[ob_var][region]), mean(all_betas[ob_var][region]), mean(high_outliers[ob_var][region]), mean(low_outliers[ob_var][region]), rank_histograms[ob_var][region].astype(float) / region_obs, region_obs, suffix, exp_name)

        pylab.clf()
        pylab.plot(times_seconds, avg_crps_values[ob_var])

        pylab.savefig("crps_avg_%s.png" % ob_var)
    return

if __name__ == "__main__":
    main()
