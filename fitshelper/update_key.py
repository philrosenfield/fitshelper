"""
Update or revert from update, fits file headers
"""
from __future__ import print_function
import argparse
import sys

from astropy.io import fits
from time import localtime, strftime


fmt = u'{0:s}: Updated {1:s} from {2:s} to {3:s}'
efmt = u'{0:s}: Doing nothing for {1:s} {2:s}'


def revert_key(key, fnames=None):
    """
    Revert an update (finds the updated key(s) in the history card)

    Parameters
    ----------
    key : string
        field name in fits header

    filenames : list of strings
        fitsfiles to change headers
    """
    fnames = fnames or []
    for fname in fnames:
        hdu = fits.open(fname, mode='update')
        hdr = hdu[0].header
        updates = [i for i in hdr['history'] if 'Updated' in i]
        update = [i for i in updates if key.upper() in i.upper()]

        if len(update) == 0:
            print('{}: no record of update to {}'.format(fname, key))
            continue

        if len(update) > 1:
            for upd in update:
                oldv, curv = map(str.strip, upd.split('from')[1].split('to'))
                if oldv != curv:
                    oldval = oldv
                    break
                else:
                    print('non-update (duplicate vals) {} {}'.format(oldv,
                                                                     curv))
        else:
            oldval = update[0].split('from')[1].split('to')[0].strip()

        update_key(key, oldval, fnames=[fname])

    return

def update_key(key, newval, fnames=None):
    """
    Change header key of fitsfile(s)

    usage: update_key.py key newval filename(s)

    Parameters
    ----------
    key : string
        field name in fits header

    newval : string
        new field name for fits header

    filenames : list of strings
        fitsfiles to change headers
    """
    fnames = fnames or []
    try:
        newval = int(newval)
    except ValueError:
        pass

    for fname in fnames:
        try:
            hdu = fits.open(fname, mode='update')
            hdr = hdu[0].header
            now = strftime("%Y-%m-%d %H:%M:%S", localtime())
            oldval = hdr[key]
            if newval != oldval:
                msg = fmt.format(now, key, oldval, newval)
                hdr[key] = newval
                hdr['history'] = msg
                print(msg)
                hdu.flush()
            else:
                print(efmt.format(fname, oldval, newval))
        except:
            print('problem with {}'.format(fname))
            print(sys.exc_info()[1])
    return


def parse_args(argv=None):
    """argparse caller for main"""
    parser = argparse.ArgumentParser(description=update_key.__doc__)

    parser.add_argument('key', type=str, help='column header to update')

    parser.add_argument('newval', type=str, help='new value')

    parser.add_argument('fnames', type=str, nargs='*',
                        help='fits file name(s)')

    return parser.parse_args(argv)


def main(argv=None):
    """main function for update_key"""
    args = parse_args(argv)

    update_key(args.key, args.newval, args.fnames)
    return

if __name__ == "__main__":
    sys.exit(main())
