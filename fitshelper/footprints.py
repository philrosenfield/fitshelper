import numpy as np
from collections import OrderedDict
from shapely.geometry import Polygon

from .utils import replace_all

def find_footprint(fitsfile, hdrext=1, closed=True):
    """
    wrapper to astropy.wcs.calc_footprint with option to append first coord
    to close the footprint.
    """
    from astropy.io import fits
    from astropy import wcs
    print(fitsfile)
    hdu = fits.open(fitsfile)
    try:
        w = wcs.WCS(hdu[hdrext].header)
    except AssertionError:
        hdrext += 1
        w = wcs.WCS(hdu[hdrext].header)
    foot = w.calc_footprint(hdu[hdrext].header)
    if closed:
        foot = np.vstack((foot, foot[0]))
    return foot

def merge_polygons(polygondic):
    pdict = {k: p for k, p in polygondic.items() if p.is_valid is True}
    first = list(pdict.keys())[0]
    ply = pdict[first]
    for i in pdict.keys():
        if i == first:
            continue
        ply = ply.union(pdict[i])
    ply2 = ply.convex_hull
    return ply2, list(pdict.keys())


def split_polygons(polygonlist, tol=49.):
    """
    Return an OrderedDict of polygons that intersect with the first value and
    an OrderedDict of polygons that do not interesect with the first value.
    """
    ins = OrderedDict()
    outs = OrderedDict()
    if len(polygonlist) > 0:
        first = list(polygonlist.keys())[0]
        ply0 = polygonlist[first].convex_hull
        ins = {first: ply0}
        for i in polygonlist.keys():
            if i == first:
                continue
            ply = polygonlist[i].convex_hull
            if ply0.intersects(ply):
                olap = ply0.intersection(ply).area / ply.area * 100
                if olap > tol:
                    ins[i] = ply
                else:
                    # print(olap)
                    outs[i] = ply
            else:
                outs[i] = ply
    return ins, outs


def group_polygons(polylist):
    """
    group a list of Polygons into ones that intersect with eachother
    returns a list of dictionaries with keys maintaining order of list index
    and values being the associated polygon
    """
    npolys = len(polylist)
    # outs = polylist
    outs = {i: p for i, p in enumerate(polylist)}
    groups = []
    while len(outs) != 0:
        ins, outs = split_polygons(outs)
        if len(ins) > 0:
            groups.append(ins)
    assert npolys == np.sum([len(g) for g in groups]), 'lost polygons'
    return groups


def parse_poly(line, closed=True, return_string=False):
    """
    Convert a polygon into N,2 np array. If closed, repeat first coords at end.
    """
    repd = {'j2000 ': '', 'gsc1 ': '', 'icrs ': '', 'multi': '',
            'polygon': '', ')': '', '(': '', 'other': ''}

    line = replace_all(line.lower(), repd).split('#')[0]
    try:
        # ds9-like format all values separated by ,
        polyline = np.array(line.strip().split(','), dtype=float)
    except:
        try:
            # shapely-like format .to_wkt(): (x0 y0, x1 y1, ...)
            xline = line.strip().replace(' ', ',') \
                                .replace(',,,', ',') \
                                .replace(',,',',')
            polyline = np.array(xline.strip().split(','), dtype=float)
        except:
            print("Do not know how to parse coords")
            sys.exit(1)

    if closed:
        if False in polyline[:2] == polyline[-2:]:
            polyline = np.append(polyline, polyline[:2])

    retv = polyline.reshape(len(polyline) // 2, 2)

    if return_string:
        retv = ','.join(['[{:.6f}, {:.6f}]'.format(*c) for c in retv])

    return retv


def parse_footprint(fname):
    """
    parse a ds9 linear footprint into a ; separated line:
    filename, polygon, central ra dec
    """
    fmt = '{};{};{}'
    with open(fname) as inp:
        lines = inp.readlines()
    # filename = fname.replace('_footprint_ds9_linear.reg', '.fits')
    polygons = [p.strip().split('#')[0] for p in lines if 'polygon' in p]
    points = [p.split('#')[0].strip() for p in lines if 'point' in p]
    texts = [p.replace('text', 'point').split('#')[0].strip()
             for p in lines if 'text' in p]
    coords = points
    if len(polygons) != len(coords):
        if len(polygons) != len(texts):
            # rolling our own...
            pcoords = [parse_poly(p) for p in polygons]
            polys = np.array([Polygon(c) for c in pcoords])
            pfmt = 'point {:.6f} {:.6f}'
            coords = [pfmt.format(p.centroid.x, p.centroid.y) for p in polys]
        else:
            coords = texts

    assert len(coords) == len(polygons), 'mismatch'

    polystr = [fmt.format(fname, polygons[i], coords[i])
               for i in range(len(coords))]
    return polystr
