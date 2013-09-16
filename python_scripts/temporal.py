
from datetime import datetime, timedelta
from itertools import chain

def _toList(iter, tolist):
    if tolist:
        return list(iter)
    else:
        return iter

class Temporal(object):
    def __init__(self, base_time, t_ens_start, t_ens_end, t_ens_step):
        self._base_time = base_time
        self._t_ens_start = t_ens_start
        self._t_ens_end = t_ens_end
        self._t_ens_step = t_ens_step
        return

    def __iter__(self):
        for t_ens in xrange(self._t_ens_start, self._t_ens_end + self._t_ens_step, self._t_ens_step):
            yield t_ens

    def __getitem__(self, index):
        return list(self)[index]

    def getTimes(self):
        return list(self)

    def getDatetimes(self, aslist=False):
        dts = ( self._base_time + timedelta(seconds=float(t_ens)) for t_ens in self )
        return _toList(dts, aslist)

    def getEpochs(self, aslist=False):
        base_epoch = (self._base_time - datetime(1970, 1, 1, 0, 0, 0)).total_seconds()
        epochs = ( base_epoch + t_ens for t_ens in self )
        return _toList(epochs, aslist)

    def getStrings(self, format, aslist=False):
        strings = ( dt.strftime(format) for dt in self.getDatetimes() )
        return _toList(strings, aslist)

class PatchedTemporal(object):
    def __init__(self, *args):
        if not all(isinstance(a, Temporal) for a in args):
            raise TypeError("All arguments to PatchedTemporal must be of type Temporal.")

        self._temporals = args
        return

    def __iter__(self):
        return chain(*self._temporals)

    def __getitem__(self, index):
        return list(self)[index]

    def getTimes(self):
        return list(self)

    def getDatetimes(self, aslist=False):
        dts = chain(*[ temp.getDatetimes() for temp in self._temporals ])
        return _toList(dts, aslist)

    def getEpochs(self, aslist=False):
        epochs = chain(*[ temp.getEpochs() for temp in self._temporals ])
        return _toList(epochs, aslist)

    def getStrings(self, format, aslist=False):
        strings = chain(*[ temp.getStrings(format) for temp in self._temporals])
        return _toList(strings, aslist)

_goshen_base_time = datetime(2009, 6, 5, 18, 0, 0)

goshen_1km_temporal = lambda start=10800, end=18000, step=300: \
    Temporal(_goshen_base_time, max(start, 10800), min(end, 18000), step)

goshen_3km_temporal = lambda start=0, end=18000, step1=1800, step2=300: \
    PatchedTemporal(
        Temporal(_goshen_base_time, max(0, start), min(10800 - step1, end), step1),
        goshen_1km_temporal(start=start, end=end, step=step2)
    )

if __name__ == "__main__":
    temp_1km = goshen_1km_temporal()
    temp_3km = goshen_3km_temporal()

    print "1km temporal testing ..."
    print type(temp_1km)
    print list(temp_1km)
    print temp_1km.getDatetimes(aslist=True)
    print temp_1km.getEpochs(aslist=True)
    for t_ens in temp_1km.getStrings("%H%M"):
        print t_ens

    print "3km temporal testing ..."
    print type(temp_3km)
    print list(temp_3km)
    print temp_3km.getDatetimes(aslist=True)
    print temp_3km.getEpochs(aslist=True)
    for t_ens in temp_3km.getStrings("%H%M"):
        print t_ens

    print "Done!"
