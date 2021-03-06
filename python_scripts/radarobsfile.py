
import numpy as np

from datetime import datetime, timedelta

from binfile import BinFile

class RadarObsFile(BinFile):
    def __init__(self, file_name, mode='r', byteorder='>'):
        super(RadarObsFile, self).__init__(file_name, mode, byteorder)

        self._readHeaders()
        self._readData()
        return

    def _readHeaders(self):
        block_size = self._read('i') # block size, I think ... 28

        timestamp = self._read('i')
        year = self._read('i')
        month = self._read('i')
        day = self._read('i')

        hour = self._read('i')
        minute = self._read('i')
        second = self._read('i')

        self._timestamp = datetime(1970, 1, 1, 0, 0, 0) + timedelta(seconds=timestamp)

        self._read('i') # previous block size? ... 28
        self._read('i') # block size, I think ... 12

        self._n_tilts = self._read('i')
        self._n_gridx = self._read('i')
        self._n_gridy = self._read('i')

        self._read('i') # previous block size? ... 12
        self._read('i') # block size, I think ... 10

        self._radar_name = self._read('10s')

        self._read('i') # previous block size? ... 10
        self._read('i') # block_size, I think ... 20

        self._radar_lat = self._read('f')
        self._radar_lon = self._read('f')
        radar_x = self._read('f')
        radar_y = self._read('f')
        self._radar_alt = self._read('f')

        self._read('i') # previous block size? ... 20
        self._read('i') # bock size, i think ... 12

        d_azimuth = self._read('f')
        range_min = self._read('f')
        range_max = self._read('f')

        self._read('i')
        return

    def _readData(self):
        self._read('i')
        elevations = self._read('f' * self._n_tilts)

        self._read('i')
        self._read('i')

        self.heights = np.transpose(self._readGrid('f', (self._n_gridx, self._n_gridy, self._n_tilts)), (2, 1, 0))

        self._read('i')
        self._read('i')

        self.range = np.transpose(self._readGrid('f', (self._n_gridx, self._n_gridy, self._n_tilts)), (2, 1, 0))

        self._read('i')
        self._variables = {}

        self._read('i')

        self._variables['vr'] = np.transpose(self._readGrid('f', (self._n_gridx, self._n_gridy, self._n_tilts)), (2, 1, 0))

        self._read('i')
        self._read('i')

        self._variables['Z'] = np.transpose(self._readGrid('f', (self._n_gridx, self._n_gridy, self._n_tilts)), (2, 1, 0))

        self._read('i')
        return

    def __getitem__(self, key):
        try:
            return self._variables[key]
        except KeyError:
            raise ValueError("Variable '%s' not found in file." % key)

if __name__ == "__main__":
#   rof = RadarObsFile("qc/manual/1km/KCYS.20090605.215309")
    rof = RadarObsFile("qc/manual/1km/KFTG.20090605.215129")
