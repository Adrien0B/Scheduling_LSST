"""Coordinate transformation methods using PyEphem
http://rhodesmill.org/pyephem/
"""

#import numpy as np
import ephem

def HESS():
    observer = ephem.Observer()
    observer.lat = "-33:27:00"
    observer.lon = "-70:40:00"
    observer.date = "2017/03/12 23:00:00"
    return observer

def transform(in_position, in_system, out_system, observer=None):
    """Simple coordinate transforms between commonly used coordinate systems.

    position: input position as a float tuple (longitude, latitude) in degrees
    in_system, out_system: coordinate system strings; one of
        horizon:    Alt,  Az
        equatorial: RA,   DEC
        galactic:   GLON, GLAT
    observer: ephem.Observer (needed when converting to / from horizon system

    Returns: transformed input position as a float tuple (longitude, latitude) in degrees
    See http://stackoverflow.com/questions/11169523/how-to-compute-alt-az-for-given-galactic-coordinate-glon-glat-with-pyephem
    """
    # Set the default observer
    observer = observer if observer else HESS()
    # Internally use radians;
    #in_position = np.radians(in_position)

    # Transform in_position to Equatorial coordinates ra, dec:
    if in_system == 'horizon':
        ra, dec = map(float, observer.radec_of(in_position[0], in_position[1]))
    elif in_system == 'equatorial':
        ra, dec = in_position
    elif in_system == 'galactic':
        galactic = ephem.Galactic(in_position[0], in_position[1])
        equatorial = ephem.Equatorial(galactic)
        ra, dec = equatorial.ra.real, equatorial.dec.real
    else:
        raise RuntimeError('in_system = %s not supported' % in_system)

    # Here we have ra, dec in radians as floats

    # Now transform Equatorial coordinates to out_system:
    if out_system == 'horizon':
        equatorial = ephem.Equatorial(ra, dec)
        body = ephem.FixedBody()
        body._ra = equatorial.ra
        body._dec = equatorial.dec
        body._epoch = equatorial.epoch
        body.compute(observer)
        out_position = body.az, body.alt
    elif out_system == 'equatorial':
        out_position = ra, dec
    elif out_system == 'galactic':
        equatorial = ephem.Equatorial(ra, dec)
        galactic = ephem.Galactic(equatorial)
        out_position = galactic.lon.real, galactic.lat.real
    else:
        raise RuntimeError('out_system = %s not supported' % out_system)

    # Clip longitude to 0 .. 360 deg range
    if out_position[0] > 360:
        out_position[0] = out_position[0] - 360

    # Return out position in radians
    return out_position
