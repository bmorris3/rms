# Licensed under the MIT License - see LICENSE

__all__ = ['Planet']

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
    "lam": 106.0,
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


class Planet(object):
    """
    Exoplanet transit properties.
    """
    def __init__(self, per=None, inc=None, a=None, t0=None, rp=None, ecc=None,
                 w=None, lam=None, b=None, duration=None, u=[0.2, 0.1],
                 limb_dark='quadratic', **kwargs):
        """
        Parameters
        ----------
        per : float
            Orbital period [days]
        inc : float
            Orbital inclination [degrees]
        a : float
            Semimajor axis [stellar radii] (a/R_star)
        t0 : float
            Mid-transit time [JD]
        rp : float
            Ratio of planet to stellar radius (Rp/Rstar)
        ecc : float
            Eccentricity
        w : float
            Argument of periastron [degrees]
        lam : float
            Projected spin-orbit angle [degrees]. ``lam=0`` is aligned,
            ``lam=90`` is perfectly misaligned.
        b : float
            Impact parameter
        duration : float
            Duration of transit [days]
        u : float
            Limb darkening parameters
        limb_dark : float
            Limb darkening formula
        """
        self.per = per
        self.inc = inc
        self.a = a
        self.t0 = t0
        self.u = u
        self.limb_dark = limb_dark
        self.rp = rp
        self.ecc = ecc
        self.w = w
        self.lam = lam
        self.b = b
        self.duration = duration

    @classmethod
    def from_hat11(cls):
        return cls(**hat11_params)
