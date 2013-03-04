
from scipy import weave
import numpy as np

def testInterp():
    interp_code = open("interp.c", 'r').read()
    pylib_code = open("pylib.c", 'r').read()

    mydata = np.arange(3 * 4 * 5, dtype=np.float32).reshape((3, 4, 5))
    z_axis = np.arange(3, dtype=np.float32).reshape((3, 1, 1)).repeat(4, axis=1).repeat(5, axis=2)
    y_axis = np.arange(4, dtype=np.float32)
    x_axis = np.arange(5, dtype=np.float32)

    slice_height = 1.51

    pts_heights = np.array([0.51, 1.51, 1.94], dtype=np.float32)
    pts_y = np.array([0.22, 0.46, 1.52], dtype=np.float32)
    pts_x = np.array([1.33, 2.89, 3.44], dtype=np.float32)

    interp_slice = np.empty((4, 5), dtype=np.float32)
    interp_pts = np.empty(pts_heights.shape, dtype=np.float32)

    interp_code = "#define INTERP_DEBUG\n%s" % interp_code

    weave.inline("interpz(mydata_array, z_axis_array, slice_height, interp_slice_array, false);", 
        ['mydata', 'z_axis', 'slice_height', 'interp_slice'], 
        support_code=pylib_code + interp_code, force=1)

    weave.inline("interppts(mydata_array, z_axis_array, pts_heights_array, y_axis_array, pts_y_array, x_axis_array, pts_x_array, interp_pts_array, false);", 
        ['mydata', 'z_axis', 'pts_heights', 'y_axis', 'pts_y', 'x_axis', 'pts_x', 'interp_pts'], 
        support_code=pylib_code + interp_code, force=1)
    print interp_slice
    print interp_pts
    return

def testCovariance():
    covariance_code = open("cov.c", 'r').read()
    pylib_code = open("pylib.c", 'r').read()

    mydata = np.arange(6 * 5 * 4 * 3).astype(np.float32).reshape((6, 5, 4, 3))
    ens_obs = np.arange(6).astype(np.float32)

    ens_cov = np.empty((5, 4, 3), dtype=np.float32)

    weave.inline("ens_covariance(mydata_array, ens_obs_array, ens_cov_array, false);",
        ['mydata', 'ens_obs', 'ens_cov'],
        support_code=pylib_code + covariance_code, 
        force=1)

    print ens_cov
    return

def testNan():
    number = 2
    result = -1

    result = weave.inline("return_val = nan(\"\");", ['number'])

    print np.isnan(result)

def main():
#   testInterp()
    testCovariance()
#   testNan()
    return

if __name__ == "__main__":
    main()
