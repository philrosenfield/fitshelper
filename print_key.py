#!/astro/apps6/anaconda2.0/bin/python
from astropy.io import fits
import sys

key = sys.argv[1]
fnames = sys.argv[2:]
for f in fnames:
    try:
        hdr = fits.getheader(f)
        print('{} {} {}'.format(f, key, hdr[key]))
    except:
        print(f, sys.exc_info()[1])
