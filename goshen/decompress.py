def decompressVariable(hd_var, dindex=None):
    if not hasattr(hd_var, 'min') or not hasattr(hd_var, 'max'):
        if dindex is None:
            return hd_var[:]
        else:
            return hd_var[dindex]

    max_i16 = np.iinfo(np.int16).max
    min_i16 = np.iinfo(np.int16).min

    if dindex is None:
        fraction = (hd_var[:].astype(np.float32) - min_i16) / (max_i16 - min_i16)

        new_shape = [ hd_var.shape[0] ]
        for idx in range(len(hd_var.shape) - 1): new_shape.append(1)
        new_shape = tuple(new_shape)

        return hd_var.min.reshape(new_shape) * (1 - fraction) + hd_var.max.reshape(new_shape) * (fraction)
    else:
        fraction = (hd_var[dindex].astype(np.float32) - min_i16) / (max_i16 - min_i16)

        if type(dindex[0]) == int:
            new_shape = [ 1 ]
        else:
            start = dindex[0].start
            if dindex[0].start is None: start = 0

            stop = dindex[0].stop 
            if dindex[0].stop is None: stop = hd_var.shape[0]

            step = dindex[0].step 
            if dindex[0].step is None: step = 1

            new_shape = [ len(range(start, stop, step)) ]

        for idx in range(len(fraction.shape) - 1): new_shape.append(1)
        new_shape = tuple(new_shape)

        decompressed_data = hd_var.min[dindex[0]].reshape(new_shape) * (1 - fraction) + hd_var.max[dindex[0]].reshape(new_shape) * (fraction)
        return decompressed_data

