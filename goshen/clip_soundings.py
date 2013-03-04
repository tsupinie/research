
import cPickle

import numpy as np

def main():
    soundings = cPickle.load(open("soundings.pkl"))

    qc = {
        'WY/NSSL2': lambda x: np.where(x > 180)[0],
        '88 and stegalllRd WY/NCAR1': lambda x: np.where(x > 490)[0],
        'NCAR2/20090605': lambda x: np.where(x > 570)[0],
        'Bushnell/NSSL1': lambda x: np.where(x < 1000)[0],
    }

    clip_sndgs = []

    for snd in soundings:
        rls = snd['release_site']

        good_idxs = qc[rls](snd['pressure'])
        for key, value in snd.iteritems():
            if type(value) == np.ndarray:
                snd[key] = value[good_idxs]
                print "Clipping %s %s ..." % (rls, key)

        clip_sndgs.append(snd)

    cPickle.dump(clip_sndgs, open("soundings_clip.pkl", 'w'), -1)
    return

if __name__ == "__main__":
    main()
