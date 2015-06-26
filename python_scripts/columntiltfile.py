
import numpy as np

from datetime import datetime, timedelta
from collections import OrderedDict

from binfile import BinFile

class ColumnTiltFile(BinFile):
    """
    ColumnTiltFile
    Author:     Tim Supinie (tsupinie@ou.edu)
    Purpose:    Read an ARPS EnKF column-tilt formatted radar observation file.
    Usage:
        ct_file = ColumnTiltFile('/path/to/file/KTLX.20110524.20000') # Read in the data
        refl = ct_file['Z'] # Get reflectivity
        rvel = ct_file['vr'] # Get radial velocity
    """

    _ctf_radar = OrderedDict([('radar_id', '4s')])
    _ctf_radar_meta1 = OrderedDict([('time_ref', 'i'), ('timestamp', 'i'), ('vcp', 'i'), ('source', 'i'), ('__dummy__', '3i'), ('nx', 'i'), ('ny', 'i'), ('nelev', 'i')])
    _ctf_run = OrderedDict([('run_name', '80s')])
    _ctf_run_meta = OrderedDict([('radar_fmt', 'i'), ('stretch_opt', 'i'), ('map_projection', 'i'), ('min_range', 'f'), ('max_range', 'f'), ('elev_type', 'i'), ('n_columns', 'i'), ('nelev', 'i'), ('__dummy__', '2i') ])
    _ctf_dom = OrderedDict([('dx', 'f'), ('dy', 'f'), ('dz', 'f'), ('dzmin', 'f'), ('center_lat', 'f'), ('center_lon', 'f'), ('std_lat1', 'f'), ('std_lat2', 'f'), ('std_lon', 'f'), ('scale_factor', 'f'),
                                  ('radar_lat', 'f'), ('radar_lon', 'f'), ('radar_alt', 'f'), ('min_elevation', 'f'), ('max_elevation', 'f')])
    _ctf_radar_over = OrderedDict([('nradover', 'i'), ('radover', '10i')])
    _ctf_dom_minmax = OrderedDict([('xmin', 'f'), ('xmax', 'f'), ('ymin', 'f'), ('ymax', 'f')])
    _ctf_time = OrderedDict([('year', 'i'), ('month', 'i'), ('day', 'i'), ('hour', 'i'), ('minute', 'i'), ('second', 'i')])
    _ctf_radar_meta2 = OrderedDict([('radar_x', 'f'), ('radar_y', 'f'), ('dazim', 'f')])

    _ctf_column_header = OrderedDict([('col_idx', 'i'), ('col_jdy', 'i'), ('col_x', 'f'), ('col_y', 'f'), ('col_lat', 'f'), ('col_lon', 'f'), ('col_sfc_hght', 'f'), ('nelev', 'i')])

    def __init__(self, file_name, mode='r', byteorder='>'):
        super(ColumnTiltFile, self).__init__(file_name, mode, byteorder)

        self._missing = -999.

        self._readHeaders()
        self._readData()
        return

    def _readHeaders(self):
        # Note: The ColumnTilt file class is aliased to _CTF at the bottom
        for block in [ _CTF._ctf_radar,
                       _CTF._ctf_radar_meta1,
                       _CTF._ctf_run,
                       _CTF._ctf_run_meta,
                       _CTF._ctf_dom,
                       _CTF._ctf_radar_over,
                       _CTF._ctf_dom_minmax,
                       _CTF._ctf_time,
                       _CTF._ctf_radar_meta2 ]:
            self._readBlock(block, self.__dict__, fortran=True)

        _elev_block = OrderedDict([('elevations', "%df" % self.nelev)])
        _elev_time_block = OrderedDict([('elev_times', "%di" % self.nelev)])

        self._readBlock(_elev_block, self.__dict__, fortran=True)
        self._readBlock(_elev_time_block, self.__dict__, fortran=True)
        return

    def _readData(self):
        columns = []
        dtype = "%df" % self.nelev
        self.vars = ['z', 'r', 'vr', 'Z'] # Wish these were encoded in the file somehow ...

        dshape = (self.nelev, self.ny, self.nx)
        self._variables = {}
        for var in self.vars:
            self._variables[var] = self._missing * np.ones(dshape, dtype=np.float32)

        for idx in range(self.n_columns):
            col_header = {}
            col_data = {}

            self._readBlock(_CTF._ctf_column_header, col_header, fortran=True)

            idx, jdy = col_header['col_idx'] - 1, col_header['col_jdy'] - 1

            for var in self.vars:
                self._readBlock(OrderedDict([(var, dtype)]), col_data, fortran=True)
                self._variables[var][:, jdy, idx] = np.array(col_data[var])

            columns.append(col_header)
        return

    def getElevationAngles(self):
        return self.elevations

    def __getitem__(self, key):
        try:
            return self._variables[key]
        except KeyError:
            raise ValueError("Variable '%s' not found in file." % key)

# Alias to _CTF to allow easier access to the format dictionaries
_CTF = ColumnTiltFile

if __name__ == "__main__":
#   rof = RadarObsFile("qc/manual/1km/KCYS.20090605.215309")
    ctf = ColumnTiltFile("/data6/tsupinie/20110524/KGLD.20110524.212909")
    print ctf['Z'].max()
    print ctf['vr'].max()
