# Licensed under the MIT License - see LICENSE
import numpy as np
import astropy.units as u

__all__ = ['Star', 'Spot']


class Star(object):
    """
    A ``batman.TransitParams``-like object for stellar parameters, to use as
    inputs for STSP.
    """
    def __init__(self, rotation_period, inc_stellar, spot_contrast, u=[0.2, 0.1],
                 limb_dark='quadratic'):
        """
        Parameters
        ----------
        rotation_period : float
            Stellar rotation period in days
        inc_stellar : float
            Stellar inclination (measured away from observer's line-of-sight)
            in units of degrees
        spot_contrast : float
            Relative intensity of a spot to the photosphere (0==perfectly dark,
            1==same as photosphere)
        u : list
            Limb darkening parameters
        limb_dark : str
            Limb darkening law
        """
        self.per = 100
        self.inc = 0
        self.a = 100
        self.t0 = 0
        self.u = u
        self.limb_dark = limb_dark
        self.rp = 0.1
        self.ecc = 0
        self.w = 90
        self.inc_stellar = inc_stellar
        self.lam = 0
        self.b = 2
        self.duration = 0.1
        self.per_rot = rotation_period
        self.spot_contrast = spot_contrast


class Spot(object):
    """
    Starspot object.
    """
    @u.quantity_input(latitude=u.deg, longitude=u.deg)
    def __init__(self, latitude, longitude, radius):
        """
        Parameters
        ----------
        latitude : `~astropy.units.Quantity`
            Spot latitude
        longitude : `~astropy.units.Quantity`
            Spot longitude
        radius : float
            Spot radius in units of stellar radii
        """
        self.latitude = latitude
        self.longitude = longitude
        self.radius = radius

    @property
    def r(self):
        return self.radius

    @property
    def theta(self):
        return np.pi/2 - self.latitude.to(u.rad).value

    @property
    def phi(self):
        return self.longitude.to(u.rad).value
