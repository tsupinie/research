
from dataIO import DataIO
import numpy as np

df = DataIO("3kmgoshen.hdf021600.01", mode="rw", format='hdf')
df.set_file_attribute('hour', 18)
df.set_file_attribute('time', 0)

df.close(out_file_name="3kmgoshenens.hdf000000")
