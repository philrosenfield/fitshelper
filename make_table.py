#!/astro/apps6/anaconda2.0/bin/python
"""
Make a table from the headers of fits files
"""
import argparse
from astropy.io import fits
import sys
import os
import numpy as np


def make_table(filelist, outfile='table.dat', clobber=False, delimiter=' ',
               nodata='...'):
    """
    Make a table from the headers of HST fits images.
    Columns are hard coded and only tested for ACS, WFC3, and WFPC2.

    Parameters
    ----------
    filelist : list
        list of fits files to read headers from
    outfile : str
        name of file to write to [table.dat]
    clobber : bool
        overwrite outfile if existing [False]
    delimiter : str
        delimiter between data [space]
        if comma, will add a line to the header to be read by MAST Portal.
    nodata : str
        fill data if the column is not found in the fits file

    Returns
    -------
        Writes or appends to outfile

    TODO
    Could be more flexible for column names.
    The formatting and the column names are separated.
    """
    head = delimiter.join('instrument detector propid target ra dec filter1 '
                          'filter2 pr_inv exptime date-obs rootname filename'
                          .split()) + '\n'

    if delimiter == ',':
        # values like exptime get string, not int, because nodata is a string.
        header = '#@string,string,string,string,ra,dec,string,string,string,string,string,string,string\n'
        header += head
        outfile = outfile.replace('dat', 'csv')
    else:
        header = '# ' + head

    acs_keys = ['INSTRUME', 'DETECTOR', 'PROPOSID', 'TARGNAME', 'RA_TARG',
                'DEC_TARG', 'FILTER1', 'FILTER2', 'PR_INV_L', 'EXPTIME',
                'DATE-OBS', 'ROOTNAME', 'EXPSTART', 'EXPEND', 'EXPTIME']

    wfc3_keys = ['INSTRUME', 'DETECTOR', 'PROPOSID', 'TARGNAME', 'RA_TARG',
                 'DEC_TARG', 'FILTER', 'FILTER2', 'PR_INV_L', 'EXPTIME',
                 'DATE-OBS', 'ROOTNAME', 'EXPSTART', 'EXPEND', 'EXPTIME']

    wfpc2_keys = ['INSTRUME', 'INSTRUME', 'PROPOSID', 'TARGNAME', 'RA_TARG',
                  'DEC_TARG', 'FILTNAM1', 'FILTNAM2', 'PROPOSID', 'EXPTIME',
                  'DATE-OBS', 'ROOTNAME', 'EXPSTART', 'EXPEND', 'EXPTIME']

    line = ''
    for filename in filelist:
        if not filename.endswith('fits'):
            print('format not supported: {}'.format(filename))
            continue
        try:
            data = fits.getheader(filename)
            inst = data['INSTRUME']
            if inst == 'WFPC2':
                keys = wfpc2_keys
            elif inst == 'ACS':
                keys = acs_keys
            else:
                keys = wfc3_keys

            fmt_list = np.array(('%({0})s %({1})s %({2})s %({3})s %({4}).6f '
                                 '%({5}).6f %({6})s %({7})s %({8})s %({9})i '
                                 '%({10})s %({11})s %({12}).6f %({13}).6f '
                                 '%({14}).6f').format(*keys).split())
            # create space for no data
            nidx = [i for i, k in enumerate(keys) if k not in data]
            # it's resonable to not have filter1, filter2, but more than that?
            if len(nidx) > 2:
                print('bad header format {}'.format(filename))

            fmt_list[nidx] = nodata
            fmt = delimiter.join(fmt_list)
        except:
            print('no header {}'.format(filename))
            continue

        line += fmt % data

        line += '{}{}\n'.format(delimiter, filename)

    if clobber or not os.path.isfile(outfile):
        wflag = 'w'
    else:
        print('appending to {}'.format(outfile))
        wflag = 'a'
        # header should already be in the file
        header = ''

    with open(outfile, wflag) as out:
        out.write(header)
        out.write(line)
        print('wrote {}'.format(outfile))


def main(argv):
    parser = argparse.ArgumentParser(
        description="Create a table from fits file header information\ne.g., make_table.py --delimiter=',' */ACS/*/*fits")

    parser.add_argument('name', nargs='*', type=str,
                        help='fits files for culling header data')

    parser.add_argument('-d', '--delimiter', type=str, default=' ',
                        help='delimiter, set to comma for MAST readable csv')

    parser.add_argument('-o', '--outfile', type=str, default='table.dat',
                        help='specify output file name, will append if output file already exists (unless -f)')

    parser.add_argument('-f', '--clobber', action='store_true',
                        help='overwite outfile if it already exists')

    args = parser.parse_args(argv)

    make_table(args.name, outfile=args.outfile, clobber=args.clobber,
               delimiter=args.delimiter)


if __name__ == '__main__':
    main(sys.argv[1:])
