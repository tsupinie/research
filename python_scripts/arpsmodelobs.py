
from binfile import BinFile

from collections import OrderedDict

import numpy as np

class ARPSModelObsFile(BinFile):
    """
    ARPSModelObsFile
    Author:     Tim Supinie (tsupinie@ou.edu)
    Purpose:    Read an ARPS EnKF model observation file (such as that created by arpsenkf or postinnov).
    Usage:
        amo_file = ARPSModelObsFile('/path/to/file/KTLXan010800', mpi_config=(nx, ny)) # Read in the data
        # The keyword argument mpi_config gives the domain decomposition configuration (nproc_x and nproc_y in arps.input). 
        # Omit the argument there is no domain decomposition.

        refl = amo_file['Z'] # Get reflectivity
        rvel = amo_file['vr'] # Get radial velocity
    """

    _amof_time = OrderedDict([('timestamp', 'i'), ('year', 'i'), ('month', 'i'), ('day', 'i'), ('hour', 'i'), ('minute', 'i'), ('second', 'i')])
    _amof_dimensions = OrderedDict([('n_tilts', 'i'), ('n_gridx', 'i'), ('n_gridy', 'i')])
    _amof_radar_name = OrderedDict([('radar_id', '4s')])
    _amof_location = OrderedDict([('radar_lat', 'f'), ('radar_lon', 'f'), ('radar_x', 'f'), ('radar_y', 'f'), ('radar_alt', 'f')])
    _amof_radar_meta = OrderedDict([('d_azim', 'f'), ('range_min', 'f'), ('range_max', 'f')])

    def __init__(self, file_name, variables=['vr', 'Z'], mpi_config=(1, 1), mode='r', byteorder='>'):
        super(ARPSModelObsFile, self).__init__(file_name, mode, byteorder)

        self._mpi = mpi_config
        self._var_names = variables
        if 'r' in mode:
            self._readHeaders()
            self._readData()
        return

    def _readHeaders(self):
        for block in [ _AMOF._amof_time,
                       _AMOF._amof_dimensions,
                       _AMOF._amof_radar_name,
                       _AMOF._amof_location,
                       _AMOF._amof_radar_meta ]:
            self._readBlock(block, self.__dict__, fortran=True)
        return

    def _readData(self):
        self.elevations = self._readGrid('f', (self.n_tilts,), fortran=True)

        block_size = self._peek('i') 

        if block_size != self.n_tilts * self.n_gridx * self.n_gridy * 4:
            self.sd_nx = self.n_gridx
            self.sd_ny = self.n_gridy

            self.n_gridx, self.n_gridy = self._subdomain2Full(self._mpi)

        self.height = self._readGrid('f', (self.n_gridx, self.n_gridy, self.n_tilts), fortran=True).T
        self.range  = self._readGrid('f', (self.n_gridx, self.n_gridy, self.n_tilts), fortran=True).T
        vars        = self._readGrid('f', (self.n_gridx, self.n_gridy, self.n_tilts, len(self._var_names)), fortran=True).T

        self._variables = dict(zip(self._var_names, vars))
        return

    def _subdomain2Full(self, mpi_config):
        nproc_x, nproc_y = mpi_config
        nx = nproc_x * (self.sd_nx - 3) + 3
        ny = nproc_y * (self.sd_ny - 3) + 3
        return nx, ny

    def __getitem__(self, var_name):
        if var_name == 'z':
            var = self.height
        elif var_name == 'r':
            var = self.range
        elif var_name in self._variables:
            var = self._variables[var_name]
        else:
            raise ValueError("Variable name %s was not given" % var_name)

        return var

# Alias to _AMOF to make accessing the format dictionaries easier.
_AMOF = ARPSModelObsFile

if __name__ == "__main__":
    ad = ARPSModelObsFile("/caps2/tsupinie/05June2009/1km-sfc-diff/KCYSan014400", mpi_config=(2, 12))
    print ad['Z'].max(), ad['Z'].min()
