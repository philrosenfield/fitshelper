"""Print fits header keys"""
import argparse
import sys

from astropy.io import fits

def print_key(key, fnames):
    """
    Print a header field of fitsfiles

    Parameters
    ----------
    key : string
        field name in fits header

    filenames : list of strings
        fits files to access
    """
    for f in fnames:
        try:
            hdr = fits.getheader(f)
            print('{0:s} {1:s} {2!s}'.format(f, key, hdr[key]))
        except:
            print(f, sys.exc_info()[1])


def parse_args(argv=None):
    """argparse caller for main"""
    parser = argparse.ArgumentParser(description=print_key.__doc__)

    parser.add_argument('key', type=str, help='field in header to print')

    parser.add_argument('fnames', type=str, nargs='*',
                        help='fits file name(s)')

    return parser.parse_args(argv)


def main(argv=None):
    """main function for update_key"""
    args = parse_args(argv)

    print_key(args.key, args.fnames)
    return

if __name__ == "__main__":
    sys.exit(main())
