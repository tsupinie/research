
import Nio as nio

import glob
from datetime import datetime, timedelta
from math import floor, ceil, fabs

def interpolate(grid, x, y, z):
    lx = int(floor(x))
    ux = int(ceil(x))
    ly = int(floor(y))
    uy = int(ceil(y))

    if ux == lx:
        if lx > 0: lx -= 1
        else:      ux += 1
    if uy == ly:
        if ly > 0: ly -= 1
        else:      uy += 1

    return ((x - lx) * (y - ly) * grid[z, ux, uy] + (ux - x) * (y - ly) * grid[z, lx, uy] +
        (x - lx) * (uy - y) * grid[z, ux, ly] + (ux - x) * (uy - y) * grid[z, lx, ly])

def main():
    files = sorted(glob.glob("3kmgoshen.hdf[0123456789]?????"))

    trajectory_points = [ [ (1, 0, 0) ],
                          [ (10, 0, 0) ],
                          [ (20, 0, 0) ],
                          [ (30, 0, 0) ],
                          [ (40, 0, 0) ] ]

    grid_spacing = 1000
    time_spacing = 300

    for file_name in files:
        hdf = nio.open_file(file_name, mode='r', format='hdf')

        u_grid = hdf.variables['u'][:]
        v_grid = hdf.variables['v'][:]

        hdf.close()

        for trajectory in trajectory_points:
            last_z, last_x, last_y = trajectory[-1]
            point_u = interpolate(u_grid, last_x, last_y, last_z)
            point_v = interpolate(v_grid, last_x, last_y, last_z)

            new_x = last_x + time_spacing * point_u / grid_spacing
            new_y = last_y + time_spacing * point_v / grid_spacing

            if new_x > u_grid.shape[1] - 1 or new_x < 0 or new_y > u_grid.shape[2] - 1 or new_y < 0:
                print "Parcel out of bounds ..."
            else:
                trajectory.append((last_z, new_x, new_y))

    print trajectory_points
    return

if __name__ == "__main__":
    main()
