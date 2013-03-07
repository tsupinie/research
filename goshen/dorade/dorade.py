
import struct
from collections import OrderedDict
import copy

import numpy as np

from binfile import BinFile
from compress import compressHRDpy, decompressHRDpy

def _error(message):
    print "Error: %s" % message
    import sys
    sys.exit()

def _warning(message):
    print "Warning: %s" % message

class DoradeFile(BinFile):
    _sswb_dict = OrderedDict([
        ('_sswb_length', 'i'), ('last_used', 'i'), ('start_time', 'i'), ('stop_time', 'i'), ('file_size', 'i'),
        ('compression_flag', 'i'), ('volume_time_stamp', 'i'), ('num_params', 'i'), ('radar_name', '8s'),
        ('start_time_float', 'd'), ('stop_time_float', 'd'), ('version_number', 'i'), ('num_key_tables', 'i'),
        ('status', 'i'), ('placeholders', 'i' * 7), ('key_table_packed', 'iii' * 8)
    ])

    _vold_dict = OrderedDict([
        ('_vold_length', 'i'), ('revision_number', 'h'), ('volume_number', 'h'), ('max_record_length', 'i'), ('project_name', '20s'),
        ('data_year', 'h'), ('data_month', 'h'), ('data_day', 'h'), ('data_hour', 'h'), ('data_minute', 'h'), ('data_second', 'h'),
        ('flight_number', '8s'), ('record_src_id', '8s'), 
        ('record_year', 'h'), ('record_month', 'h'), ('record_day', 'h'),
        ('n_sens_descr', 'h')
    ])

    _radd_dict = OrderedDict([
        ('_radd_length', 'i'), ('name', '8s'), ('radar_constant', 'f'), ('peak_power', 'f'), ('noise_power', 'f'), ('receiver_gain', 'f'), ('antenna_gain', 'f'), ('system_gain', 'f'), ('horiz_beam_width', 'f'), ('vert_beam_width', 'f'),
        ('radar_type', 'h'), ('scan_mode', 'h'), ('antenna_rot_vel', 'f'), ('scan_param1', 'f'), ('scan_param2', 'f'),
        ('n_param_descr', 'h'), ('n_additional_descr', 'h'),
        ('data_compression', 'h'), ('data_reduction', 'h'), ('data_reduction_param1', 'f'), ('data_reduction_param2', 'f'),
        ('radar_longitude', 'f'), ('radar_latitude', 'f'), ('radar_altitude', 'f'),
        ('unambig_velocity', 'f'), ('unambig_range', 'f'), 
        ('n_frequencies', 'h'), ('n_interpulse_per', 'h'), ('frequencies', 'f' * 5), ('interpulse_per', 'f' * 5),

        # Extended block
        ('extension', 'i'), ('config_name', '8s'), ('config_num', 'i'), 
        ('aperature_size', 'f'), ('field_of_view', 'f'), ('aperature_eff', 'f'), 
        ('ext_frequencies', 'f' * 11), ('ext_interpulse_per', 'f' * 11), ('pulse_width', 'f'),
        ('primary_baseline', 'f'), ('secondary_baseline', 'f'), ('pc_tx_bandwidth', 'f'), ('pc_type', 'i'),
        ('site_name', '20s')
    ])

    _lidr_dict = OrderedDict([
        ('_lidr_length', 'i'), ('lidar_name', '8s'), ('lidar_constant', 'f'), ('pulse_energy', 'f'), ('peak_power', 'f'), ('pulse_width', 'f'), ('aperature_size', 'f'), ('field_of_view', 'f'), ('aperature_eff', 'f'), ('beam_divergence', 'f'),
        ('lidar_type', 'h'), ('scan_mode', 'h'), ('mirror_rot_vel', 'f'), ('lidar_scan_param1', 'f'), ('lidar_scan_param2', 'f'),
        ('n_param_desc', 'h'), ('n_total_descr', 'h'),
        ('data_compression', 'h'), ('data_reduction', 'h'), ('data_reduction_param1', 'f'), ('data_reduction_param2', 'f'),
        ('lidar_longitude', 'f'), ('lidar_latitude', 'f'), ('lidar_altitude', 'f'),
        ('unambig_velocity', 'f'), ('unambig_range', 'f'), 
        ('n_wavelengths', 'i'), ('pulse_rep_freq', 'f'), ('frequencies', 'f' * 10),
    ])

    _parm_dict = OrderedDict([
        ('_parm_length', 'i'), ('name', '8s'), ('description', '40s'), ('units', '8s'), ('interpulse_used', 'h'), ('frequency_used', 'h'),
        ('receiver_bandwidth', 'f'), ('pulse_width', 'h'), ('polarization', 'h'), ('n_samples', 'h'), ('binary_format', 'h'),
        ('threshold_param', '8s'), ('threshold_value', 'f'), 
        ('scale_factor', 'f'), ('bias_factor', 'f'), ('bad_data_flag', 'i'),

        # Extended block
        ('extension', 'i'), ('config_name', '8s'), ('config_num', 'i'), ('offset_to_data', 'i'), ('mks_conversion', 'f'), 
        ('n_qnames', 'i'), ('qnames', '32s'), ('n_criteria', 'i'), ('crit_names', '32s'), 
        ('n_cells', 'i'), ('meters_to_first_cell', 'f'), ('meters_between_cells', 'f'), ('eff_unambig_velocity', 'f')
    ])

    _csfd_dict = OrderedDict([
        ('_csfd_length', 'i'), ('n_segments', 'i'), ('meters_to_first', 'f'), ('meters_between', 'f' * 8), ('n_cells', 'h' * 8)
    ])

    _cfac_dict = OrderedDict([
        ('_cfac_length', 'i'), ('azimuth', 'f'), ('elevation', 'f'), ('range_delay', 'f'), ('longitude', 'f'), ('latitude', 'f'), ('pressure_alt', 'f'), ('physical_alt', 'f'), 
        ('platform_u', 'f'), ('platform_v', 'f'), ('platform_w', 'f'), ('platform_heading', 'f'), ('platform_roll', 'f'), ('platform_pitch', 'f'), ('platform_drift', 'f'), 
        ('rotation_angle', 'f'), ('tilt_angle', 'f')
    ])

    _swib_dict = OrderedDict([
        ('_swib_length', 'i'), ('sweep_comment', '8s'), ('sweep_number', 'i'), ('n_rays', 'i'), ('true_start_angle', 'f'), ('true_end_angle', 'f'), ('fixed_angle', 'f'), ('filter_flag', 'i')
    ])

    _ryib_dict = OrderedDict([
        ('_ryib_length', 'i'), ('sweep_number', 'i'), ('julian_day', 'i'), ('hour', 'h'), ('minute', 'h'), ('second', 'h'), ('millisecond', 'h'), ('azimuth', 'f'), ('elevation', 'f'), 
        ('peak_tx_power', 'f'), ('scan_rate', 'f'), ('ray_status', 'i')
    ])

    _asib_dict = OrderedDict([
        ('_asib_length', 'i'), ('longitude', 'f'), ('latitude', 'f'), ('altitude_msl', 'f'), ('altitude_agl', 'f'),
        ('antenna_u', 'f'), ('antenna_v', 'f'), ('antenna_w', 'f'), ('heading', 'f'), ('roll', 'f'), ('pitch', 'f'), ('drift', 'f'),
        ('beam_sweep_angle', 'f'), ('beam_scan_angle', 'f'),
        ('air_u', 'f'), ('air_v', 'f'), ('air_w', 'f'), ('heading_chg_rate', 'f'), ('pitch_chg_rate', 'f')
    ])

    _rktb_dict = OrderedDict([
        ('_rktb_length', 'i'), ('angle2index', 'i'), ('n_indices', 'i'), ('ray_lut_offset', 'i'), ('angle_table_offset', 'i'), ('n_rays', 'i')
    ])

    def __init__(self, file_name, mode, byteorder='<'):
        super(DoradeFile, self).__init__(file_name, mode, byteorder)
        self._mode = mode
        self._defaults = []
        self._defaults = dir(self)

        if 'r' in mode:
            self._readHeaders()
            self._readSweep()
            self._readBackMatter()
        return

    def _readHeaders(self):
        self._readCommentBlock()
        self._readSuperSweepIDBlock()
        self._readVolumeDescriptor()
        self._readSensorDescriptors()
        return

    def _readCommentBlock(self):
        _marker = self._peek('4s')
        if _marker != "COMM":
            return

        _marker = self._read('4s')
        self._comm_size = self._read('i')
        self.comment = self._read("%ds" % (self._comm_size - 8))
        return

    def _readSuperSweepIDBlock(self):
        _marker = self._read('4s')
        if _marker != "SSWB": 
            _error("Expected volume descriptor marker (SSWB), but got '%s'" % _marker)

        block_length = self._peek('i')
        block_mod = OrderedDict([ (key, value) for key, value in DoradeFile._sswb_dict.iteritems() ])
        if block_length == 200:
            block_mod['radar_name'] = '12s'

        self._readBlock(block_mod, self.__dict__)
        self.key_table = self._unpackTable(self.key_table_packed, 3)

#       print self.placeholders
#       print self.key_table
        return

    def _readVolumeDescriptor(self):
        _marker = self._read('4s')
        if _marker != "VOLD": 
            _error("Expected volume descriptor marker (VOLD), but got '%s'" % _marker)

        self._readBlock(DoradeFile._vold_dict, self.__dict__)

        return

    def _readSensorDescriptors(self):
        self.radar_descriptors = []
#       for sens_descr in range(self.n_sens_descr):
        while self._peek('4s') == "RADD":
            descriptor = {}

            _marker = self._read('4s')
            if _marker != "RADD": 
                _error("Expected sensor descriptor marker (RADD), but got '%s'" % _marker)

            block_length = self._peek('i')
            block_cutoff = len(DoradeFile._radd_dict)
            if block_length == 144:
                block_cutoff = 30

            block_dict = OrderedDict([ (key, value) for idx, (key, value) in enumerate(DoradeFile._radd_dict.items()) if idx < block_cutoff ])
            self._readBlock(block_dict, descriptor)


            lidr_descriptor = {}
            _marker = self._peek('4s')
            if _marker == "LIDR":
                _marker = self._read('4s')
                self._readBlock(DoradeFile._lidr_dict, lidr_descriptor)

            # Figure out what to do with this lidar block ...

            descriptor['parameters']  = self._readParameterDescriptors(descriptor['n_param_descr'])
            descriptor['cell_range']  = self._readCellRangeInfo()
            descriptor['corr_factor'] = self._readCorrectionFactorDescriptor()

            self.radar_descriptors.append(descriptor)

        if len(self.radar_descriptors) != self.n_sens_descr:
            _warning("Number of sensor descriptors does not match the number in the file.  Fixing this.")
            self.n_sens_descr = len(self.radar_descriptors)
        return

    def _readParameterDescriptors(self, n_parameter_descr):
        parameter_descriptors = []

        for parm_desc in range(n_parameter_descr):
            descriptor = {}

            _marker = self._read('4s')
            if _marker != "PARM": 
                _error("Expected parameter descriptor marker (PARM), but got '%s'" % _marker)


            block_length = self._peek('i')
            block_cutoff = len(DoradeFile._parm_dict)
            if block_length == 104:
                block_cutoff = 16

            block_dict = OrderedDict([ (key, value) for idx, (key, value) in enumerate(DoradeFile._parm_dict.items()) if idx < block_cutoff ])
            self._readBlock(block_dict, descriptor)

            parameter_descriptors.append(descriptor)

        return parameter_descriptors

    def _readCellRangeInfo(self):
        descriptor = {}
        _marker = self._read('4s')
        if _marker != "CELV" and _marker != "CSFD": 
            _error("Expected cell range vector marker (CELV) or cell spacing table marker (CSFD), but got '%s'" % _marker)

        if _marker == "CELV":
            descriptor['_celv_length']  = self._read('i')
            descriptor['vector_length'] = self._read('i')
            descriptor['vector']        = self._read('f' * ((descriptor['_celv_length'] - 12) / 4))
        elif _marker == "CSFD":
            self._readBlock(DoradeFile._csfd_dict, descriptor)

        return descriptor

    def _readCorrectionFactorDescriptor(self):
        descriptor = {}
        _marker = self._peek('4s')
        if _marker != "CFAC":
            return
#           _error("Expected correction factor descriptor marker (CFAC), but got '%s'" % _marker)

        _marker = self._read('4s')
        self._readBlock(DoradeFile._cfac_dict, descriptor)
        return descriptor

    def _readSweep(self):
        self._readSweepDescriptor()
        self._readRays(self.radar_descriptors[0]['parameters'])
        return

    def _readSweepDescriptor(self):
        _marker = self._read('4s')
        if _marker != "SWIB": 
            _error("Expected sweep descriptor marker (SWIB), but got '%s'" % _marker)

        self._readBlock(DoradeFile._swib_dict, self.__dict__)
        return

    def _readRays(self, parameter_descriptors):
        self._rays = []
        from datetime import datetime, timedelta

        read_time = timedelta(0)
        decompress_time = timedelta(0)

        for ray in range(self.n_rays):
            descriptor = {}
            descriptor['ray_info']      = self._readRayInfoBlock()

            platform_info = self._readPlatformInfoBlock()
            if platform_info != {}:
                descriptor['platform_info'] = platform_info

            descriptor['param_data']    = {}
            for param_desc in parameter_descriptors:
                _marker = self._read('4s')

                if _marker != "RDAT":
                    _error("Expected radar data marker (RDAT) for ray %d, but got '%s'" % (ray, _marker))

                radar_data_length = self._read('i')
                parameter_name = self._read('8s')

                if parameter_name != param_desc['name']:
                    _error("Expected parameter %s, but got %s" % (param_desc['name'], parameter_name))

                data_type = {1:'b', 2:'h', 3:'i', 4:'f'}[ param_desc['binary_format'] ]
                data_width = {1:1, 2:2, 3:4, 4:4}[ param_desc['binary_format'] ]

                type_str = data_type * ((radar_data_length - 16) / data_width)
                begin = datetime.now()
                data_compressed = np.array(self._read(type_str))
                read_time += (datetime.now() - begin)

                begin = datetime.now()
                data = decompressHRDpy(data_compressed)
                decompress_time += (datetime.now() - begin)

                data_remapped = self._remapFromFile(data, param_desc)
                descriptor['param_data'][parameter_name] = data_remapped

            self._rays.append(descriptor)
#       print "Time to read data:", read_time
#       print "Time to decompress data:", decompress_time
        return

    def _remapFromFile(self, data, parameter_desc):
        return np.where(data > -10000, data / parameter_desc['scale_factor'] - parameter_desc['bias_factor'], data)

    def _readRayInfoBlock(self):
        descriptor = {}
        _marker = self._read('4s')
        if _marker != "RYIB": 
            _error("Expected ray info block marker (RYIB), but got '%s'" % _marker)

        self._readBlock(DoradeFile._ryib_dict, descriptor)
        return descriptor

    def _readPlatformInfoBlock(self):
        descriptor = {}

        _marker = self._peek('4s')
        if _marker != "ASIB":
            return descriptor
#           _error("Expected platform info block marker (ASIB), but got '%s'" % _marker)

        _marker = self._read('4s')
        self._readBlock(DoradeFile._asib_dict, descriptor)
        return descriptor

    def _readBackMatter(self):
        _marker = self._read('4s')
        if _marker != "NULL":
            _error("Expected null block marker (NULL), but got '%s'" % _marker)
        null_size = self._read('i')

        _marker = self._read('4s')
        if _marker != "RKTB":
            _error("Expected rotation angle table block marker (RKTB), but got '%s'" % _marker)

        self._readBlock(DoradeFile._rktb_dict, self.__dict__)

#       for rktb_key in DoradeFile._rktb_dict.keys():
#           print rktb_key, self.__dict__[rktb_key]

        rot_table_length = self._rktb_length - struct.calcsize("".join(DoradeFile._rktb_dict.values())) - struct.calcsize('i') * (self.n_indices + 1)
        rot_table_entry_length = struct.calcsize('fii')

        self._rot_lut = self._read('i' * (self.n_indices + 1))
        self._rot_table = self._unpackTable(self._read('fii' * (rot_table_length / rot_table_entry_length)), 3)

#       print self._rot_lut
#       print self._rot_table

        _marker = self._peek('4s')
        if _marker == "SEDS":
            _marker = self._read('4s')
            self._seds_length = self._read('i')
            self.editor_history = self._read("s%d" % (self._seds_length - 8))
#           _error("Expected editor history block marker (SEDS), but got '%s'" % _marker)

        return

    def close(self):
        if 'w' in self._mode:
            self._writeHeaders()
            self._writeSweep()
            self._writeBackMatter()

            # Need to update key table offsets here (the first gets self._rktb_offset and the second gets self._seds_offset)
            self.key_table[0][0] = self._rktb_offset

            if hasattr(self, '_seds_offset'):
                self.key_table[1][0] = self._seds_offset

            offset = struct.calcsize("".join(DoradeFile._sswb_dict.values()[:-1])) + 4 # +4 for the marker
            self._seek(offset)
            self._write(self._packTable(self.key_table), '24i')

        super(DoradeFile, self).close()

    def _writeHeaders(self):
        self._writeCommentBlock()
        self._writeSuperSweepIDBlock()
        self._writeVolumeDescriptor()
        self._writeSensorDescriptors()
        return

    def _writeCommentBlock(self):
        if hasattr(self, 'comment'):
            self._write("COMM", '4s')
            self._write(len(self.comment) + 8, 'i')
            self._write(self.comment, "%ds" % len(self.comment))
        return

    def _writeSuperSweepIDBlock(self):
        self.key_table_packed = self._packTable(self.key_table)

        block_mod = OrderedDict([ (key, value) for key, value in DoradeFile._sswb_dict.iteritems() ])
        if self._sswb_length == 200:
            block_mod['radar_name'] = '12s'

        self._write("SSWB", '4s')
        self._writeBlock(block_mod, self.__dict__)
        return

    def _writeVolumeDescriptor(self):
        self._write("VOLD", '4s')
        self._writeBlock(DoradeFile._vold_dict, self.__dict__)
        return

    def _writeSensorDescriptors(self):
        for n_desc in range(self.n_sens_descr):
            descriptor = self.radar_descriptors[n_desc]

            self._write("RADD", '4s')

            block_cutoff = len(DoradeFile._radd_dict)
            if descriptor['_radd_length'] == 144:
                block_cutoff = 30

            block_dict = OrderedDict([ (key, value) for idx, (key, value) in enumerate(DoradeFile._radd_dict.items()) if idx < block_cutoff ])
            self._writeBlock(block_dict, descriptor)

            self._writeParameterDescriptors(descriptor['parameters'])
            self._writeCellRangeInfo(descriptor['cell_range'])
            self._writeCorrectionFactorDescriptor(descriptor['corr_factor'])
        return

    def _writeParameterDescriptors(self, param_descrs):
        for descriptor in param_descrs:
            self._write("PARM", '4s')

            block_cutoff = len(DoradeFile._parm_dict)
            if descriptor['_parm_length'] == 104:
                block_cutoff = 16

            block_dict = OrderedDict([ (key, value) for idx, (key, value) in enumerate(DoradeFile._parm_dict.items()) if idx < block_cutoff ])
            self._writeBlock(block_dict, descriptor)
        return

    def _writeCellRangeInfo(self, cell_range_info):
        if '_celv_length' in cell_range_info:
            self._write("CELV", '4s')
            self._write(cell_range_info['_celv_length'], 'i')
            self._write(cell_range_info['vector_length'], 'i')
            self._write(cell_range_info['vector'], 'f' * ((cell_range_info['_celv_length'] - 12) / 4))
        else:
            self._write("CSFD", '4s')
            self._writeBlock(DoradeFile._csfd_dict, cell_range_info)
        return

    def _writeCorrectionFactorDescriptor(self, corr_factor):
        if hasattr(self, '_cfac_length'):
            self._write("CFAC", '4s')
            self._writeBlock(DoradeFile._cfac_dict, corr_factor)
        return

    def _writeSweep(self):
        self._writeSweepDescriptor()
        self._writeRays(self.radar_descriptors[0]['parameters'])
        return

    def _writeSweepDescriptor(self):
        self._write("SWIB", '4s')
        self._writeBlock(DoradeFile._swib_dict, self.__dict__)
        return

    def _writeRays(self, parameters):
        from datetime import datetime, timedelta
        write_time = timedelta(0)
        compress_time = timedelta(0)

        for idx, ray in enumerate(self._rays):
            self._rays[idx]['offset'] = self._tell()
            self._writeRayInfoBlock(ray['ray_info'])

            if 'platform_info' in ray:
                self._writePlatformInfoBlock(ray['platform_info'])

            for param in parameters:
                data_type = {1:'b', 2:'h', 3:'i', 4:'f'}[ param['binary_format'] ]
                data_width = {1:1, 2:2, 3:4, 4:4}[ param['binary_format'] ]

                parameter_name = param['name']
                data = ray['param_data'][parameter_name]
                data_remap = self._remapToFile(data, param)

                begin = datetime.now()
                data_compressed = compressHRDpy(data_remap)
                compress_time += (datetime.now() - begin)

                block_length = data_width * len(data_compressed) + 16

                self._write("RDAT", '4s')
                self._write(block_length, 'i')
                self._write(parameter_name, '8s')

                begin = datetime.now()
                self._write(data_compressed, data_type * len(data_compressed))
                write_time += (datetime.now() - begin)

            self._rays[idx]['length'] = self._tell() - self._rays[idx]['offset']

#       print "Time to write data:", write_time
#       print "Time to compress data:", compress_time
        return


    def _remapToFile(self, data, parameter_desc):
        return np.where(data > -10000, (data + parameter_desc['bias_factor']) * parameter_desc['scale_factor'], data)

    def _writeRayInfoBlock(self, ray_info_block):
        self._write("RYIB", '4s')
        self._writeBlock(DoradeFile._ryib_dict, ray_info_block)
        return

    def _writePlatformInfoBlock(self, platform_info_block):
        if platform_info_block != {}:
            self._write("ASIB", '4s')
            self._writeBlock(DoradeFile._asib_dict, platform_info_block)
        return

    def _writeBackMatter(self):
        self._write("NULL", '4s')
        self._write(8, 'i')

        self._rktb_offset = self._tell()
        rktb_header_size = struct.calcsize("".join(DoradeFile._rktb_dict.values()))

        rot_lut_length = self.n_indices + 1
#       self.ray_lut_offset = self._rktb_offset + rktb_header_size + 4 # +4 for the marker
#       self.angle_table_offset = self.ray_lut_offset + rot_lut_length

        self._write("RKTB", '4s')
        self._writeBlock(DoradeFile._rktb_dict, self.__dict__)

        # Need to update the rotation table for offsets and lengths of the rays
        for row in range(len(self._rays)):
            self._rot_table[row][1] = self._rays[row]['offset']
            self._rot_table[row][2] = self._rays[row]['length']

        rot_table_length = self._rktb_length - struct.calcsize("".join(DoradeFile._rktb_dict.values())) - struct.calcsize('i') * rot_lut_length
        rot_table_entry_length = struct.calcsize('fii')

        self._write(self._rot_lut, 'i' * rot_lut_length)
        self._write(self._packTable(self._rot_table), 'fii' * (rot_table_length / rot_table_entry_length))

        if hasattr(self, 'editor_history'):
            self._seds_offset = self._tell()
            self._write("SEDS", '4s')
            self._write(len(self.editor_history) + 8, 'i')
            self._write(self.editor_history, "s%d" % len(self.editor_history))
        return

    def _packTable(self, table):
        table_packed = [ cell for row in table for cell in row ]
        return table_packed

    def _unpackTable(self, table_packed, n_cols):
        n_rows = len(table_packed) / n_cols
        table = [ [ table_packed[nr * n_cols + nc] for nc in range(n_cols) ] for nr in range(n_rows) ]
        return table

    def copyValues(self, other):
        for attr in dir(other):
            if attr not in other._defaults:
                setattr(self, attr, copy.copy(getattr(other, attr)))
        return

    def getSweep(self, parameter_name):
        param_descriptor = [ p for p in self.radar_descriptors[0]['parameters'] if p['name'] == parameter_name ][0]
        threshold_param = param_descriptor['threshold_param']
        threshold_value = param_descriptor['threshold_value']

        sweep = np.empty((self.n_rays, len(self._rays[0]['param_data'][parameter_name])))
        for idx in range(self.n_rays):
            if threshold_param != "":
                threshold_data = self._rays[idx]['param_data'][threshold_param.strip()]
                sweep[idx] = np.where(threshold_data < threshold_value, self._rays[idx]['param_data'][parameter_name], -int("0x8000", 0))
            else:
                sweep[idx] = self._rays[idx]['param_data'][parameter_name]

        return sweep

    def setSweep(self, parameter_name, sweep):
        if sweep.shape[0] != self.n_rays:
            _warning("Number of rays in the sweep is different than in the file.  Nothing done.")
            return

        for idx in range(self.n_rays):
            self._rays[idx]['param_data'][parameter_name] = sweep[idx]
        return

    def getCellRanges(self):
        return np.array(self.radar_descriptors[0]['cell_range']['vector'])

    def getRayAngles(self):
        end_angle = self.true_end_angle
        if end_angle < self.true_start_angle:
            end_angle += 360
        angles = np.linspace(self.true_start_angle, end_angle, self.n_rays)
        return np.where(angles >= 360, angles - 360, angles)

def open_file(file_name, mode, byteorder='<'):
    return DoradeFile(file_name, mode, byteorder)

if __name__ == "__main__":
    from datetime import datetime

    begin = datetime.now()
    dor = open_file("/data6/tsupinie/goshen/raw/KCYS/swp.1090605213925.KCYS.275.0.5_SUR_v409", 'r')#, byteorder='>')
    print "Time to read file:", datetime.now() - begin
    dor.close()

    dor2 = open_file("swp.copy", 'w')
    dor2.copyValues(dor)
    begin = datetime.now()
    dor2.close()
    print "Time to write file:", datetime.now() - begin

    dor3 = open_file("swp.copy", 'r')
    dor3.close()

    import sys; sys.exit()

    cell_range = dor3.getCellRanges()
    start_range_idx = np.argmin(np.abs(0 - cell_range))
    end_range_idx = np.argmin(np.abs(1000000 - cell_range))

    sweep = dor3.getSweep("DZ")
    xs, ys = np.meshgrid(np.arange(sweep.shape[0]), np.arange(sweep.shape[1]))
    import pylab
    pylab.pcolormesh(xs, ys, sweep.T)
    pylab.colorbar()

    pylab.show()

    import glob
#   files = glob.glob("raw/MWR05XP/*/swp.1090605215744.MWR_05XP*")
#   files.sort(key=lambda x: float(x.split("/")[2]))
#   for file in files:
#       dor = open_file(file, 'r')
#       print file
#       print "Horizontal Beam Width (HBW): ", dor.radar_descriptors[0]['horiz_beam_width']
#       print "Vertical Beam Width (VBW):   ", dor.radar_descriptors[0]['vert_beam_width']
#       print "VBW / HBW:                   ", dor.radar_descriptors[0]['vert_beam_width'] / dor.radar_descriptors[0]['horiz_beam_width']
#
#       print "Start bearing:               ", dor.true_start_angle
#       print "End bearing:                 ", dor.true_end_angle
#       print "Number of rays:              ", dor.num_rays
#       print "Ray spacing:                 ", (360 + dor.true_start_angle - dor.true_end_angle) / (dor.num_rays - 1)
#       print
