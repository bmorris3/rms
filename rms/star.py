# Licensed under the MIT License - see LICENSE

from .planet import Planet

__all__ = ['Star']


class Star(object):
    """
    An object for stellar parameters, to use as inputs for STSP.
    """
    def __init__(self, planet=None, rotation_period=None, inc_stellar=None,
                 spot_contrast=0.7, u=[0.2, 0.1], rho_s=1.0):
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
        u : float
            Quadratic limb darkening parameters
        limb_dark : float
            Limb darkening formula
        planet : `~rms.Planet`
            Planet parameters. If planet is None, a non-transiting planet will
            be used for STSP computations.
        rho_s : float
            Stellar density in units of the solar density
        """
        self.inc_stellar = inc_stellar
        self.per_rot = rotation_period
        self.spot_contrast = spot_contrast

        if planet is None:
            planet = Planet.non_transiting()

        self.planet = planet
        self.u = u
        self.rho_s = rho_s
