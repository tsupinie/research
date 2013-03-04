import numpy as np

import struct
from operator import mul

class BinFile(object):
    def __init__(self, file_name, mode='r', byteorder="<"):
        self._bin_file = open(file_name, "%sb" % mode)
        self._byteorder = byteorder
        self._file_pos = self._bin_file.tell()

        self._read_target = False
        return

    def close(self):
        self._bin_file.close()
        self._file_pos = -1
        return

    def _read(self, type_string, _peeking=False):
        self._bin_file.seek(self._file_pos)

        size = struct.calcsize(type_string)
        raw = self._bin_file.read(size)

        if raw == "":
            return ""

        data = struct.unpack("%s%s" % (self._byteorder, type_string), raw)

        if not _peeking:
            self._file_pos = self._bin_file.tell()

        if len(data) == 1:
            if type_string[-1] == 's':
                return data[0].strip("\0")
            else:
                return data[0]
        else:
            return list(data)

    def _write(self, value, type_string):
        if type(value) not in [ list, tuple ]:
            value = [ value ]

        self._bin_file.write(struct.pack("%s%s" % (self._byteorder, type_string), *value))
        return

    def _readGrid(self, type_string, shape):
        length = reduce(mul, shape)
        return np.array(self._read(type_string * length)).reshape(shape, order='F')

    def _readBlock(self, type_dict, dest_dict):
        for key, type in type_dict.iteritems():
            dest_dict[key] = self._read(type)
        return

    def _writeBlock(self, type_dict, src_dict):
        for key, type in type_dict.iteritems():
            self._write(src_dict[key], type)
        return

    def _peek(self, type_string):
        return self._read(type_string, _peeking=True)

    def _seek(self, location, anchor=0):
        self._bin_file.seek(location, anchor)
        return

    def _tell(self):
        return self._bin_file.tell()

    def _ateof(self):
        return self._bin_file.read(1) == ""

if __name__ == "__main__":
    bf = BinFile("qc/manual/1km/KCYS.20090605.215744")
    assert bf._peek('i') == 28
    assert bf._read('i') == 28

    print "Unit tests done."
