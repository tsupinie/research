
from binfile import BinFile

import numpy as np

class ARPSModelObsFile(BinFile):
    def __init__(self, file_name, variables=['vr', 'Z'], mpi_config=(1, 1), mode='r', byteorder='>'):
        super(ARPSModelObsFile, self).__init__(file_name, mode, byteorder)

        self._mpi = mpi_config
        self._var_names = variables
        if 'r' in mode:
            self._readHeaders()
            self._readData()
        return

    def _readHeaders(self):
        block_size = self._read('i')
        self._read('i' * (block_size / 4)) # Supposed to be time, worthless as far as I can tell
        assert block_size == self._read('i')

        block_size = self._read('i')
        self._ntilts, sd_nx, sd_ny = self._read('i' * (block_size / 4))
        nproc_x, nproc_y = self._mpi
        self._nx = nproc_x * (sd_nx - 3) + 3
        self._ny = nproc_y * (sd_ny - 3) + 3

        assert block_size == self._read('i')

        block_size = self._read('i')
        self._radarid = self._read("%ds" % block_size)
        assert block_size == self._read('i')

        block_size = self._read('i')
        self._radar_lat, self._radar_lon, self._radar_x, self._radar_y, self._radar_alt = self._read('f' * (block_size / 4))
        assert block_size == self._read('i')

        block_size = self._read('i')
        self._dazim, self._range_min, self._range_max = self._read('f' * (block_size / 4))
        assert block_size == self._read('i')

        block_size = self._read('i')
        self._elev_angles = self._read('f' * (block_size / 4))
        assert block_size == self._read('i')
        return

    def _readData(self):
        block_size = self._read('i') 
        assert block_size == self._ntilts * self._nx * self._ny * 4
        self._height = np.transpose(self._readGrid('f', (self._nx, self._ny, self._ntilts)), [2, 1, 0])
        assert block_size ==  self._read('i')

        block_size = self._read('i')
        assert block_size == self._ntilts * self._nx * self._ny * 4
        self._range = np.transpose(self._readGrid('f', (self._nx, self._ny, self._ntilts)), [2, 1, 0])
        assert block_size == self._read('i')

        self._variables = np.transpose(self._readGrid('f', (self._nx, self._ny, self._ntilts, len(self._var_names))), [3, 2, 1, 0])

        return

    def __getitem__(self, var_name):
        try:
            var_idx = self._var_names.index(var_name)
        except ValueError:
            raise ValueError("Variable name %s was not given" % var_name)

        return self._variables[var_idx]

if __name__ == "__main__":
    ad = ARPSModelObsFile("/caps1/tsupinie/3km-fixed-radar/KCYSan014400", (2, 12))
    print ad['Z'].max()
