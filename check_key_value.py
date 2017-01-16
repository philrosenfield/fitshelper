#!/astro/apps6/anaconda2.0/bin/python
from astropy.io import fits
import sys

key = sys.argv[1]
value = sys.argv[2]
fnames = sys.argv[3:]
print('These files have the following {}={}'.format(key, value))
for f in fnames:
    try:
        hdr = fits.getheader(f)
        if hdr[key].lower() == value.lower():
            print('{}'.format(f))
    except:
        print(sys.exc_info()[1])
