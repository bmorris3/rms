# Licensed under the MIT License - see LICENSE

__all__ = ['Star']


class Star(object):
    """
    A ``batman.TransitParams``-like object for stellar parameters, to use as
    inputs for STSP.
    """
    def __init__(self, planet=None, rotation_period=None, inc_stellar=None,
                 spot_contrast=0.7):
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
        planet : `~rms.Planet`
            Planet parameters
        """

        self.inc_stellar = inc_stellar
        self.per_rot = rotation_period
        self.spot_contrast = spot_contrast
        self.planet = planet
