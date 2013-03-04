
from util import goshen_3km_proj, goshen_3km_gs, drawPolitical, setupMapProjection

import numpy as np

import matplotlib
from matplotlib.patches import Circle
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap

import Nio as nio

def makeComposite(radar_data, type='max'):
    if type == "max":
        return np.array(radar_data).max(axis=0)
    else:
        return

def makePlot(radar_data, gs, map, title, file_name, radars=None):
    gs_x, gs_y = gs
    nx, ny = radar_data.shape
    xs, ys = np.meshgrid(gs_x * np.arange(nx), gs_y * np.arange(ny))

    pylab.contourf(xs, ys, radar_data, levels=np.arange(10, 80, 10))
    pylab.colorbar()

    if radars:
        for radar_id, (radar_lat, radar_lon) in radars.iteritems():
            radar_x, radar_y = map(radar_lon, radar_lat)
            pylab.gca().add_patch(Circle([radar_x, radar_y], 230000, fc='none', ec='#333333'))

    drawPolitical(map, scale_len=75)

    pylab.suptitle(title)
    pylab.savefig(file_name)
    pylab.close()
    return

def main():
    files = [ 'hdf/KCYS/3km/manual/goshen.hdfrefl2d010691', 'hdf/KRIW/3km/goshen.hdfrefl2d010785', 'hdf/KFTG/3km/goshen.hdfrefl2d010638' ]
    radar_data = [ nio.open_file(f, format='hdf', mode='r').variables['refl2d'][0] for f in files ]
    radars = { 'KCYS':(41.15194, -104.80611), 'KFTG':(39.78667, -104.54583), 'KRIW':(43.06611, -108.47722) }

    composite = makeComposite(radar_data)

    proj = setupMapProjection(goshen_3km_proj, goshen_3km_gs)
    map = Basemap(**proj)

    makePlot(composite, goshen_3km_gs, map, "Radar Obs composite (2100 UTC)", "radar_composite_010800.png", radars=radars)
#   for fn, rd in zip(files, radar_data):
#       rid = fn[4:8]
#       makePlot(rd, goshen_3km_gs, map, "Radar Obs composite component (%s, 2100 UTC)" % rid, "radar_composite_comp_%s_010800.png" % rid)
    return

if __name__ == "__main__":
    main()
