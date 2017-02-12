"""Cross match a MAST Discovery Portal table (with s_regions) to a source table"""
from __future__ import print_function
import argparse
import os
import sys

from collections import OrderedDict
import numpy as np
import pandas as pd

from fitshelper.footprints import merge_polygons, parse_poly, group_polygons
from shapely.geometry import Polygon, Point

from ..utils import fixnames


def cross_match(lit_cat, mast_cat, plot=False, ra='RAJ2000',
                dec='DEJ2000', namecol='SimbadName'):
    """
    Cross match a literature catalog (of points) with a MAST Discovery
    Portal catalog (with polygons as s_regions).

    Parameters
    ----------
    lit_cat : str
        filename of the literature catalog

    mast_cat : str
        filename of the MAST Discovery Portal catalog

    plot : bool [False]
        make a plot of the matched catalogs.
        (probably better to upload them back to MAST)

    ra : str [RAJ2000]
        ra column name

    dec : str [DEJ2000]
        dec column name

    namecol : str [SimbadName]
        name of column from lit_cat to add to new matched catalog

    Returns
    -------
    file: [mast_cat]_matched_[lit_cat].csv
    A culled mast_cat to matches to lit_cat
    will have columns "group" and "litidx"
        group : s_regions that are likely covering (>49%) the same area
        litidx : the matched index (line not counting headers) from the lit_cat.
    """
    if plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()

    # this code can also be used to cross match the HST query result
    # against a MAST catalog. So the "matched" catalog will be read in
    # and has a different header than the one downloaded directly from
    # the Discovery Portal.
    # This try/except attepts to raise a KeyError if the header isn't
    # read correctly.
    try:
        mast = pd.read_csv(mast_cat, header=4)
        mast['s_region']
    except KeyError:
        mast = pd.read_csv(mast_cat,  header=0, skiprows=[1])
        mast['s_region']

    # if group has already been added to the table, this table is being used
    # to cross match another table.
    gstr = 'group'
    if gstr in mast.columns:
        gstr += '1'
    mast[gstr] = np.nan

    lstr = 'litidx'
    if lstr in mast.columns:
        lstr += '1'
    mast[lstr] = np.nan

    try:
        lit = pd.read_csv(lit_cat,  header=0, skiprows=[0])
        lit[ra]
    except KeyError:
        lit = pd.read_csv(lit_cat,  header=0, skiprows=[0])
        lit[ra]

    lit[gstr] = np.nan

    if namecol in mast.columns:
        # the lit_cat should be the basis of the output file.
        savelit = True
        mast = fixnames(mast)
        lit[namecol] = np.nan
    else:
        savelit = False
        lit = fixnames(lit)

    radecs = [Point(r,d) for r,d in zip(lit[ra], lit[dec])]

    plys = [Polygon(parse_poly(mast['s_region'].iloc[i]))
            for i in range(len(mast['s_region']))]

    # Group observations by overlapping footprints
    grouped_plys = group_polygons(plys)

    # For the following search, make one polygon per group
    group = OrderedDict()
    for i, g in enumerate(grouped_plys):
        p, inds = merge_polygons(g)
        group[i] = ({'p': p, 'inds': inds})
        # Assign group id to mast table
        mast[gstr].iloc[inds] = i

    # loop through all combinations
    for g, pdict in group.items():
        for i, radec in enumerate(radecs):
            # Ignore single filter observations (1-item groups)
            if len(pdict['inds']) == 1 or \
                len(np.unique(mast.iloc[pdict['inds']].filters)) < 2:
                continue
            # Is the radec point in the s_region?
            if pdict['p'].contains(radec):
                if np.isfinite(lit[gstr].iloc[i]):
                    # there are complex regions that are not grouped together
                    # so more than one radec point can be in a polygon.
                    # there has to be a better way than making a single item
                    # string...
                    tmp = lit[gstr].iloc[i]
                    lit[gstr].iloc[i] = ','.join(['{:g}'.format(t)
                                                     for t in [tmp, g]])
                # add group id to literature table
                lit[gstr].iloc[i] = g
                if savelit:
                    if len(np.unique(mast.iloc[pdict['inds']][namecol])) > 1:
                        print('More than one lit source in a s_region. Not supported')
                    else:
                        lit[namecol].iloc[i] = np.unique(mast.iloc[pdict['inds']][namecol])[0]
                # add literature index to the mast table
                mast[lstr].iloc[pdict['inds']] = i

                if plot:
                    p = np.array(radec.to_wkt()
                                      .replace('POINT (', '')
                                      .replace(')','').split(), dtype=float)
                    v = parse_poly(pdict['p'].to_wkt())
                    ax.plot(v[:, 0], v[:, 1])
                    ax.plot(p[0], p[1], 'o')

    fins, = np.nonzero(np.isfinite(mast[lstr]))
    print('{} lit values, {} MAST values, {} matches.'.format(len(radecs),
                                                              len(plys),
                                                              len(fins)))
    if savelit:
        df = lit
    else:
        df = mast.loc[fins]
        if namecol in lit.columns:
            df[namecol] = lit.iloc[df[lstr]][namecol].tolist()
        else:
            print('{0:s} not found in lit cat, will not add a name column'
                  .format(namecol))

    fname = \
        '{}_matched_{}'.format(os.path.split(mast_cat.replace('.csv',''))[1],
                               os.path.split(lit_cat)[1])
    df.to_csv(fname, index=False)
    #write out lit.to_csv?

    if plot:
        plt.savefig('cross_match.png')
        plt.close()


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="cross match a literature catalog with a MAST catalog")

    parser.add_argument('-p', '--plot', action='store_true',
                        help='make a plot')

    parser.add_argument('--ra', type=str, default='RAJ2000',
                        help='literature ra column name')

    parser.add_argument('--dec', type=str, default='DEJ2000',
                        help='literature dec column name')

    parser.add_argument('lit_cat', type=str,
                        help='literature catalog')

    parser.add_argument('mast_cat', type=str,
                        help='MAST catalog (or the one with s_region)')

    parser.add_argument('--namecol', type=str, default='SimbadName',
                        help='Column from lit to add to MAST')

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    cross_match(args.lit_cat, args.mast_cat, plot=args.plot, ra=args.ra,
                dec=args.dec, namecol=args.namecol)


if __name__ == "__main__":
    sys.exit(main())
