# Licensed under the MIT License - see LICENSE

__all__ = ['Star']

hat11_params = {
    "a": 14.770858177139598,
    "fp": None,
    "b": 0.143443142704645,
    "ecc": 0.2744922585429323,
    "inc_stellar": 100.0,
    "rho_star": 1.8082748494218275,
    "limb_dark": "quadratic",
    "per": 4.8878025843894006,
    "per_rot": 29.192083459347486,
    "lam": 0.0, #106.0,
    "t_secondary": None,
    "t0": 2454605.8914623754,
    "w": 18.03890135536712,
    "duration": 0.097973065981468405,
    "rp": 0.058330305324663184,
    "u": [
        0.6407001070602456,
        0.047761746059854178
    ],
    "inc": 89.34708509841991
}


class Star(object):
    """
    A ``batman.TransitParams``-like object for stellar parameters, to use as
    inputs for STSP.
    """
    def __init__(self, rotation_period=None, inc_stellar=None,
                 spot_contrast=0.7, u=[0.2, 0.1], limb_dark='quadratic',
                 planet=False):
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
        planet : bool
            Drop in the planet parameters for HAT-P-11 b. Default is False.
        """
        if planet:
            for attr in hat11_params.keys():
                setattr(self, attr, hat11_params[attr])

        else:
            self.per = 100
            self.inc = 0
            self.a = 100
            self.t0 = 0
            self.u = u
            self.limb_dark = limb_dark
            self.rp = 0.1
            self.ecc = 0
            self.w = 90
            self.lam = 0
            self.b = 2
            self.duration = 0.1

        if inc_stellar is not None:
            self.inc_stellar = inc_stellar
        if rotation_period is not None:
            self.per_rot = rotation_period
        if spot_contrast is not None:
            self.spot_contrast = spot_contrast
