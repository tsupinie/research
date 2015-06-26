
import numpy as np

from datetime import datetime, timedelta
from collections import OrderedDict

from binfile import BinFile

class GridTiltFile(BinFile):
    """
    GridTiltFile
    Author:     Tim Supinie (tsupinie@ou.edu)
    Purpose:    Read an ARPS EnKF grid-tilt formatted radar observation file.
    Usage:
        gt_file = GridTiltFile('/path/to/file/KTLX.20110524.20000') # Read in the data
        refl = gt_file['Z'] # Get reflectivity
        rvel = gt_file['vr'] # Get radial velocity
    """
    _gtf_time = OrderedDict([('timestamp', 'i'), ('year', 'i'), ('month', 'i'), ('day', 'i'), ('hour', 'i'), ('minute', 'i'), ('second', 'i')])
    _gtf_dimensions = OrderedDict([('n_tilts', 'i'), ('n_gridx', 'i'), ('n_gridy', 'i')])
    _gtf_name = OrderedDict([('radar_name', '10s')])
    _gtf_location = OrderedDict([('radar_lat', 'f'), ('radar_lon', 'f'), ('radar_x', 'f'), ('radar_y', 'f'), ('radar_alt', 'f')])
    _gtf_radar_meta = OrderedDict([('d_azimuth', 'f'), ('range_min', 'f'), ('range_max', 'f')])

    def __init__(self, file_name, mode='r', byteorder='>'):
        super(GridTiltFile, self).__init__(file_name, mode, byteorder)

        self._readHeaders()
        self._readData()
        return

    def _readHeaders(self):
        for block in [ _GTF._gtf_time,
                       _GTF._gtf_dimensions,
                       _GTF._gtf_name,
                       _GTF._gtf_location,
                       _GTF._gtf_radar_meta ]:
            self._readBlock(block, self.__dict__, fortran=True)

        self.timestamp = datetime(1970, 1, 1, 0, 0, 0) + timedelta(seconds=self.timestamp)
        return

    def _readData(self):
        self.elevations = self._readGrid('f', (self.n_tilts,), fortran=True)
        self.heights    = self._readGrid('f', (self.n_gridx, self.n_gridy, self.n_tilts), fortran=True).T
        self.range      = self._readGrid('f', (self.n_gridx, self.n_gridy, self.n_tilts), fortran=True).T

        self._variables = {}
        self._variables['vr'] = self._readGrid('f', (self.n_gridx, self.n_gridy, self.n_tilts), fortran=True).T
        self._variables['Z']  = self._readGrid('f', (self.n_gridx, self.n_gridy, self.n_tilts), fortran=True).T
        return

    def __getitem__(self, key):
        try:
            return self._variables[key]
        except KeyError:
            raise ValueError("Variable '%s' not found in file." % key)

#Alias to _GTF to make it easier to access the format dictionaries
_GTF = GridTiltFile

if __name__ == "__main__":
    gtf = GridTiltFile("/data6/tsupinie/goshen/qc/manual/1km/KFTG.20090605.215129")
    print gtf['Z'].max(), gtf['Z'].min()
