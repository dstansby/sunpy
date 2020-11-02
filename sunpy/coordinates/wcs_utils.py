from astropy.coordinates import SkyCoord, HeliocentricMeanEcliptic
import astropy.units as u
import astropy.wcs.utils
from astropy.wcs import WCS

from sunpy.util import MetaDict

from .frames import (
    BaseCoordinateFrame,
    Heliocentric,
    HeliographicCarrington,
    HeliographicStonyhurst,
    Helioprojective,
    SunPyBaseCoordinateFrame,
)

__all__ = ['solar_wcs_frame_mapping', 'solar_frame_to_wcs_mapping',
           'observer_coordinate']


def solar_wcs_frame_mapping(wcs):
    """
    This function registers the coordinates frames to their FITS-WCS coordinate
    type values in the `astropy.wcs.utils.wcs_to_celestial_frame` registry.

    Parameters
    ----------
    wcs : astropy.wcs.WCS

    Returns
    -------
    astropy.coordinates.BaseCoordinateFrame
    """

    if hasattr(wcs, "coordinate_frame"):
        return wcs.coordinate_frame

    dateobs = wcs.wcs.dateobs or None

    # SunPy Map adds 'heliographic_observer' and 'rsun' attributes to the WCS
    # object. We check for them here, and default to None.
    if hasattr(wcs.wcs, 'aux'):
        observer = observer_coordinate(wcs)
    elif hasattr(wcs, 'heliographic_observer'):
        observer = wcs.heliographic_observer
    else:
        observer = None

    if hasattr(wcs, 'rsun'):
        rsun = wcs.rsun
    else:
        rsun = None

    # Truncate the ctype to the first four letters
    ctypes = {c[:4] for c in wcs.wcs.ctype}

    if {'HPLN', 'HPLT'} <= ctypes:
        return Helioprojective(obstime=dateobs, observer=observer, rsun=rsun)

    if {'HGLN', 'HGLT'} <= ctypes:
        return HeliographicStonyhurst(obstime=dateobs)

    if {'CRLN', 'CRLT'} <= ctypes:
        return HeliographicCarrington(obstime=dateobs, observer=observer)

    if {'SOLX', 'SOLY'} <= ctypes:
        return Heliocentric(obstime=dateobs, observer=observer)


def solar_frame_to_wcs_mapping(frame, projection='TAN'):
    """
    For a given frame, this function returns the corresponding WCS object.
    It registers the WCS coordinates types from their associated frame in the
    `astropy.wcs.utils.celestial_frame_to_wcs` registry.

    Parameters
    ----------
    frame : astropy.coordiantes.BaseCoordinateFrame
    projection : str, optional

    Returns
    -------
    astropy.wcs.WCS
    """
    wcs = WCS(naxis=2)

    if hasattr(frame, 'rsun'):
        wcs.rsun = frame.rsun
    else:
        wcs.rsun = None

    if hasattr(frame, 'observer') and isinstance(frame.observer, BaseCoordinateFrame):
        wcs.heliographic_observer = frame.observer
    else:
        wcs.heliographic_observer = None

    if isinstance(frame, SunPyBaseCoordinateFrame):

        if frame.obstime:
            wcs.wcs.dateobs = frame.obstime.utc.isot

        if isinstance(frame, Helioprojective):
            xcoord = 'HPLN' + '-' + projection
            ycoord = 'HPLT' + '-' + projection
            wcs.wcs.cunit = ['arcsec', 'arcsec']
        elif isinstance(frame, Heliocentric):
            xcoord = 'SOLX'
            ycoord = 'SOLY'
            wcs.wcs.cunit = ['deg', 'deg']
        elif isinstance(frame, HeliographicCarrington):
            xcoord = 'CRLN' + '-' + projection
            ycoord = 'CRLT' + '-' + projection
            wcs.wcs.cunit = ['deg', 'deg']
        elif isinstance(frame, HeliographicStonyhurst):
            xcoord = 'HGLN' + '-' + projection
            ycoord = 'HGLT' + '-' + projection
            wcs.wcs.cunit = ['deg', 'deg']

    else:
        return None

    wcs.wcs.ctype = [xcoord, ycoord]

    return wcs


astropy.wcs.utils.WCS_FRAME_MAPPINGS.append([solar_wcs_frame_mapping])
astropy.wcs.utils.FRAME_WCS_MAPPINGS.append([solar_frame_to_wcs_mapping])


class ObserverSystem:
    """
    An observer coordinate system supported by sunpy.

    Parameters
    ----------
    frame : astropy.coordinates.BaseFrame
        The coordinate frame.
    fits_keywords : list of str
        The list of three FITS keywords used by this system.
    coord_kwargs : list of str
        The list of keyword argument names used to construct the coordinate
        SkyCoord. These must be in the same order as ``fits_keywords``.
    units : list of `astropy.units.Unit`
        The units associated with the value stored in each keyword.
        These must be in the same order as ``fits_keywords``.
    aux_kw : list of str, optional
        List of keyword names used by the WCS aux attribute. See `astropy.wcs.Auxprm`
        for a list of allowed values.
    kwargs : dict, optional
        Any additional keyword arguments are handed to the SkyCoord construcor
        when creating the coordinate.
    """
    def __init__(self, frame, fits_keywords, coord_kwargs, units, aux_kw=None, **kwargs):
        self.frame = frame
        self.fits_keywords = fits_keywords
        self.coord_kwargs = coord_kwargs
        self.units = units

        if aux_kw is None:
            aux_kw = self.fits_keywords
        self.aux_kw = aux_kw

        self.kwargs = kwargs

    def missing_keys(self, meta):
        """
        Return a list of the FITS keywords missing in meta.
        """
        return [keyword for keyword in self.fits_keywords if keyword not in meta]

    def coord(self, meta):
        """
        Returns the observer coordinate constructed from meta.
        """
        obstime = meta.get('date-obs', None)
        args = [meta[fkw] for fkw in self.fits_keywords]
        sc = SkyCoord(*args, obstime=obstime, frame=self.frame, unit=self.units, **self.kwargs)

        # We need to specially handle an observer location provided in Carrington
        # coordinates.  To create the observer coordinate, we need to specify the
        # frame, but defining a Carrington frame normally requires specifying the
        # frame's observer.  This loop is the problem.  Instead, since the
        # Carrington frame needs only the Sun-observer distance component from the
        # frame's observer, we create the same frame using a fake observer that has
        # the same Sun-observer distance.
        if isinstance(sc.frame, HeliographicCarrington):
            fake_observer = HeliographicStonyhurst(0*u.deg, 0*u.deg, sc.radius,
                                                   obstime=sc.obstime)
            fake_frame = sc.frame.replicate(observer=fake_observer)
            hgs = fake_frame.transform_to(HeliographicStonyhurst(obstime=sc.obstime))

            # HeliographicStonyhurst doesn't need an observer, but adding the observer
            # facilitates a conversion back to HeliographicCarrington
            return SkyCoord(hgs, observer=hgs)
        return sc.transform_to('heliographic_stonyhurst')


# A mapping of coordinate systems to the metadata keys that have to be present
# for the observer coordinate to be fully defined
supported_observer_systems = [ObserverSystem(HeliocentricMeanEcliptic,
                                             ['haex_obs', 'haey_obs', 'haez_obs'],
                                             ['x', 'y', 'z'],
                                             [u.m, u.m, u.m],
                                             representation_type='cartesian'),
                              ObserverSystem(HeliographicStonyhurst,
                                             ['hgln_obs', 'hglt_obs', 'dsun_obs'],
                                             ['lon', 'lat', 'radius'],
                                             [u.deg, u.deg, u.m]),
                              ObserverSystem(HeliographicCarrington,
                                             ['crln_obs', 'crlt_obs', 'dsun_obs'],
                                             ['lon', 'lat', 'radius'],
                                             [u.deg, u.deg, u.m],
                                             aux_kw=['crln_obs', 'hglt_obs', 'dsun_obs']),
                              ]

# The same as above, but just systems supported by the WCS auxillary metadata
supported_aux_systems = supported_observer_systems[1:3]


def _match_system(meta):
    missing_meta = {}
    for system in supported_observer_systems:
        missing = system.missing_keys(meta)
        if len(missing) == 0:
            return system

        missing_meta[system.frame] = missing

    # Haven't found a fully matching system, collect missing keys and error
    error_message = "".join(
        [f"For frame '{frame}' the following metadata is missing: "
         f"{','.join(keys)}\n" for frame, keys in missing_meta.items()])
    raise ValueError(error_message)


def observer_coordinate(meta):
    """
    Return the observer coordinate defined by some FITS metadata.

    Parameters
    ----------
    meta : sunpy.util.MetaDict, astropy.wcs.WCS
    """
    if isinstance(meta, WCS):
        aux_dict = {}
        aux = meta.wcs.aux
        print(aux)
        for system in supported_aux_systems:
            for fits_kw, aux_kw in zip(system.fits_keywords, system.aux_kw):
                aux_dict[fits_kw] = getattr(aux, aux_kw)
        aux_dict['date-obs'] = meta.wcs.dateobs or None
        meta = MetaDict(aux_dict)

    system = _match_system(meta)
    return system.coord(meta)
