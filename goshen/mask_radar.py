
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

import glob
import os

from dorade.dorade import DoradeFile

def plotSweep(ranges, angles, refl):
    rs, thetas = np.meshgrid(ranges, np.radians(angles))

    xs = rs * np.cos(np.pi / 2 - thetas)
    ys = rs * np.sin(np.pi / 2 - thetas)

    print xs.shape, ys.shape, refl.shape

    print xs
    print ys

    pylab.pcolor(xs, ys, refl, vmin=10, vmax=80, cmap=matplotlib.cm.get_cmap('jet'))
#   pylab.pcolor(xs, ys, mask, vmin=0, vmax=1, cmap=matplotlib.cm.get_cmap('binary'), alpha=0.1)
    pylab.xlim([-100000, 100000])
    pylab.ylim([-100000, 100000])
    pylab.savefig("sweep.png")

def main():
    files_qc = [ glob.glob("raw/KCYS/swp.*0.5_SUR_v*")[30] ]
    files_noband = [ glob.glob("raw/KCYS_noband/swp.*0.5_SUR_v*")[30] ]
    for file_qc, file_nb in zip(files_qc, files_noband):
        assert os.path.basename(file_qc) == os.path.basename(file_nb)

        print os.path.basename(file_qc)

        dor_qc = DoradeFile(file_qc, 'r')
        dor_qc.close()

#       print dor_qc.radar_descriptors[0]['parameters']

        dor_nb = DoradeFile(file_nb, 'r', byteorder='>')
        dor_nb.close()

        bad_raw_refl = dor_qc.getSweep('REF')
        qc_refl = dor_qc.getSweep('DZ')
        good_raw_refl = dor_nb.getSweep('REF')

        n_rays, n_gates = good_raw_refl.shape
        cell_range = np.array(dor_qc.getCellRanges())[:(n_gates + 1)]
        angles = dor_qc.getRayAngles()

        good_qc_refl = np.where(bad_raw_refl == qc_refl, good_raw_refl, qc_refl)

        cutoff = np.argmin(np.abs(cell_range - 120000))
        good_qc_refl[:, :cutoff] = qc_refl[:, :cutoff]

#       plotSweep(cell_range, angles, good_qc_refl)

        dor_good = DoradeFile(os.path.basename(file_qc), 'w')
        dor_good.copyValues(dor_qc)
        dor_good.setSweep('DZ', good_qc_refl)
        dor_good.close()
    return

if __name__ == "__main__":
    main()
