
import numpy as np
from scipy import weave

def decompressHRDpy(compressed_data, debug=False):
    decompressed_data = []
    idx = 0

    if debug:
        print compressed_data

    while idx < len(compressed_data) and compressed_data[idx] != 1:

        count = compressed_data[idx] & int("0x7fff", 0)
        good_data = compressed_data[idx] & int("0x8000", 0)

        if debug:
            print count, bool(good_data)

        if good_data:
            decompressed_data.extend(compressed_data[(idx + 1):(idx + count + 1)])
            idx += count + 1
        else:
            decompressed_data.extend([ -int("0x8000", 0) for jdy in range(count) ])
            idx += 1

    return np.array(decompressed_data)

def decompressHRD(compressed_data):
    return

def compressHRDpy(decompressed_data):
    compressed_data = []
    missing = -int("0x8000", 0)
    idx = 0

    while idx < len(decompressed_data):
        run_end_idx = idx

        if idx < len(decompressed_data) - 1 and decompressed_data[idx] == missing and decompressed_data[idx + 1] == missing:
            # Run of missing data
            while run_end_idx < len(decompressed_data) and decompressed_data[run_end_idx] == missing: 
                run_end_idx += 1

            compressed_data.append(run_end_idx - idx)
            idx = run_end_idx
        else:
            # Run of good data
            while run_end_idx < len(decompressed_data) and (run_end_idx - idx < 2 or decompressed_data[run_end_idx] != missing or 
                (run_end_idx < len(decompressed_data) - 1 and decompressed_data[run_end_idx] == missing and decompressed_data[run_end_idx + 1] != missing)): 
                run_end_idx += 1
        
            compressed_data.append(missing + (run_end_idx - idx))
            compressed_data.extend(decompressed_data[idx:run_end_idx].astype(int))

            idx = run_end_idx

    compressed_data.append(1)
    return compressed_data
