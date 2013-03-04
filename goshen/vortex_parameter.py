
import Nio as nio

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab
import matplotlib.transforms as transforms
from matplotlib.scale import LogScale

import cPickle

from plot_sounding import plotSounding

def intersect(*args):
    if len(args) == 2:
        return np.intersect1d(*args)
    else:
        return np.intersect1d(args[0], intersect(*args[1:]))

def vortexParameter(u, v, normalize=True):
    new_shape = [ 4 ]
    new_shape.extend(u.shape)

    new_shape[-1] -= 2
    new_shape[-2] -= 2

    new_shape = tuple(new_shape)

    def vortexComponent(comp, normalize=True):
        if normalize:
            comp_avg = comp[1:-1, 1:-1]
        else:
            comp_avg = np.zeros(new_shape[1:])

        comp_norm = np.zeros(new_shape)
        comp_norm[0] = comp[:-2, 1:-1]
        comp_norm[1] = comp[1:-1, :-2]
        comp_norm[2] = comp[2:, 1:-1]
        comp_norm[3] = comp[1:-1, 2:]

        return comp_norm

    u_norm = vortexComponent(u, normalize)
    v_norm = vortexComponent(v, normalize)

    return u_norm[0] - v_norm[1] - u_norm[2] + v_norm[3]

def main():
    z_index = 1
    subset = (slice(80, 130), slice(100, 150))
    exp_name = "1kmnam"
    ens_name = "1kmnam"
    time = 14400
    x_snd, y_snd = 115, 90

    hdf = nio.open_file("%s/%s.hdf%06d" % (exp_name, ens_name, time), mode='r', format='hdf')

    u = hdf.variables['u'][z_index]
    v = hdf.variables['v'][z_index]
    w = hdf.variables['w'][z_index]
    u_clip = u[1:-1, 1:-1]
    v_clip = v[1:-1, 1:-1]
    w_clip = w[1:-1, 1:-1]

    u_snd = hdf.variables['u'][:, y_snd, x_snd]
    v_snd = hdf.variables['v'][:, y_snd, x_snd]
    pt_snd = hdf.variables['pt'][:, y_snd, x_snd]
    p_snd = hdf.variables['p'][:, y_snd, x_snd]
    qv_snd = hdf.variables['qv'][:, y_snd, x_snd]

    vortex_param = vortexParameter(u, v, normalize=False)

    x, y = np.meshgrid(np.arange(vortex_param.shape[1]) + 1, np.arange(vortex_param.shape[0]) + 1)
    pylab.contourf(x[subset], y[subset], vortex_param[subset])
    pylab.colorbar()
    pylab.quiver(x[subset], y[subset], u_clip[subset], v_clip[subset])
    pylab.savefig("%s.vp.%06d.png" % (ens_name, time))

    pylab.clf()
    vorticity = (np.gradient(v, 1000)[1] - np.gradient(u, 1000)[0])[1:-1, 1:-1]
    pylab.contourf(x[subset], y[subset], vorticity[subset])
    pylab.colorbar()
    pylab.quiver(x[subset], y[subset], u_clip[subset], v_clip[subset])
    pylab.savefig("%s.vort.%06d.png" % (ens_name, time))

    pylab.clf()
    divergence = (np.gradient(u, 1000)[1] + np.gradient(v, 1000)[0])[1:-1, 1:-1]
    pylab.contourf(x[subset], y[subset], divergence[subset])
    pylab.colorbar()
    pylab.quiver(x[subset], y[subset], u_clip[subset], v_clip[subset])
    pylab.savefig("%s.div.%06d.png" % (ens_name, time))

    pylab.clf()
    pylab.contourf(x[subset], y[subset], w_clip[subset])
    pylab.colorbar()
    pylab.quiver(x[subset], y[subset], u_clip[subset], v_clip[subset])
    pylab.savefig("%s.w.%06d.png" % (ens_name, time))

    plotSounding("%s.snd-%03d,%03d.%06d.png" % (ens_name, x_snd, y_snd, time), pt=pt_snd, p=(p_snd / 100.), qv=qv_snd, u=u_snd, v=v_snd)

    proximity = cPickle.load(open("proximity2.pkl", 'r'))
    u_obs_snd = proximity['u_wind']
    v_obs_snd = proximity['v_wind']
    t_obs_snd = proximity['temperature']
    td_obs_snd = proximity['dewpoint']
    p_obs_snd = proximity['pressure']

    good = (u_obs_snd != 9999.0) & (v_obs_snd != 9999.0) & (t_obs_snd != 999.0) & (td_obs_snd != 999.0) & (p_obs_snd != 9999.0)
    good_idxs = np.where(good)[0]

    plotSounding("proximity.bushnell.052155.png", p=p_obs_snd[good_idxs], t=t_obs_snd[good_idxs], td=td_obs_snd[good_idxs], u=u_obs_snd[good_idxs], v=v_obs_snd[good_idxs])
    return

if __name__ == "__main__":
    main()
