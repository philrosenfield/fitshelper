import numpy as np
from shapely.geometry import Polygon


def replace_all(text, dic):
    """perfrom text.replace(key, value) for all keys and values in dic"""
    for old, new in dic.items():
        text = text.replace(old, new)
    return text


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


def parse_poly(line, closed=True, return_string=False):
    """
    Convert a polygon into N,2 np array. If closed, repeat first coords at end.
    """
    repd = {'J2000 ': '', 'GSC1 ': '', 'ICRS ': '', 'multi': '',
            'polygon': '', ')': '', '(': ''}

    line = replace_all(line.lower(), repd).split('#')[0]
    try:
        # ds9-like format all values separated by ,
        polyline = np.array(line.strip().split(','), dtype=float)
    except:
        # shapely-like format (x0 y0, x1 y1, ...)
        line = line.strip().replace(' ', ',').replace(',,', ',')
        polyline = np.array(line.strip().split(','), dtype=float)

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
