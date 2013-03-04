
import Nio as nio

import glob

def main():
    base_path = "/caps1/tsupinie/1km-control-no-ua"
    files = glob.glob("%s/ena???.hdf0*" % base_path)

    for file in files:
        hdf = nio.open_file(file, mode='r', format='hdf')

        for var in ['u', 'v', 'w', 'pt', 'p', 'qv']:
            if var not in hdf.variables:
                print "%s incomplete ... " % file
                break

        hdf.close()
    return

if __name__ == "__main__":
    main()
