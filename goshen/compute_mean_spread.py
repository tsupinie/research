
import Nio as nio

import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

def spread(ens_members, variable):
    member_shape = ens_members[0].variables[variable].shape
    full_shape = [ len(ens_members) ]
    full_shape.extend(member_shape)
    full_shape = tuple(full_shape)

    full_data = np.empty(full_shape)
    for n_ens, hdf in enumerate(ens_members):
        full_data[n_ens] = hdf.variables[variable][:]

    return full_data.std(axis=0, ddof=1)

def plotSpread(spread, title, file_name):
    pylab.clf()
    ny, nx = spread.shape
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    pylab.contourf(x, y, spread)
    pylab.colorbar()
    pylab.title(title)
    pylab.savefig(file_name)
    return

def main():
    n_ens_members = 40
    bounds = (slice(3, -3), slice(12, -12), slice(12, -12))

    ens_members = []
    for n_ens in range(n_ens_members):
        file_name = "1km-control-20120617/ena%03d.hdf014400" % (n_ens + 1)
        hdf = nio.open_file(file_name, mode='r', format='hdf')
        ens_members.append(hdf)

    for var in ['u', 'v', 'w', 'pt', 'p', 'qv']:
        var_spread = spread(ens_members, var)[bounds]
        print "Max spread in %s: %f" % (var, var_spread.max())
        print "Min spread in %s: %f" % (var, var_spread.min())
        print "Median spread in %s: %f" % (var, np.median(var_spread))
        print "Mean spread in %s: %f" % (var, var_spread.mean())

        plotSpread(var_spread[12], "Spread in %s at k=13 at 2100 UTC 05 June 2009" % var, "sprd-%s.010800.png" % var)

    return

if __name__ == "__main__":
    main()
