# Licensed under the MIT License - see LICENSE
import numpy as np

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

hat7_params = {
    "a": 4.1502,
    "fp": None,
    "b": 0.4978,
    "ecc": 0.0,
    "inc_stellar": 0,  # Lund 2014
    "rho_star": 0.18875,  # Lund 2014
    "limb_dark": "quadratic",
    "per": 2.204737,  # Morris
    "per_rot": 13.0,  # Lund 2014
    "lam": 155,  # Albrecht et al 2012
    "t_secondary": None,
    "t0": 2454954.357463,
    "w": 90,
    "duration": 0.13701750368545595,  # Morris
    "rp": 0.077590, # Morris
    "u": [
        0.3525, # Morris
        0.168
    ],
    "inc": 83.111  # Morris
}

null_planet_params = {
    "a": 100,
    "fp": None,
    "b": 2,
    "ecc": 0.0,
    "per": 100,
    "t0": 2454954.357463,
    "w": 90,
    "rp": 0.001,
    "inc": 0,
    "lam": 0,
    "duration": 0
}


class Planet(object):
    """
    Exoplanet transit properties.
    """
    def __init__(self, per=None, inc=None, a=None, t0=None, rp=None, ecc=None,
                 w=None, lam=None, b=None, duration=None, **kwargs):
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
            Impact parameter. If none is supplied and ``ecc=0``, the impact
            parameter will be computed for you.
        duration : float
            Duration of transit [days]. If none is supplied and ``ecc=0``, one
            will be computed for you.
        """
        self.per = per
        self.a = a
        self.t0 = t0
        self.rp = rp
        self.ecc = ecc
        self.w = w
        self.lam = lam

        if inc != 0:
            # If given ecc and inc but not b, compute b
            if b is None and ecc == 0 and inc is not None:
                b = a * np.cos(np.radians(inc))

            # If given ecc and b but not inc, compute inc
            if b is not None and ecc == 0 and inc is None:
                inc = np.degrees(np.arccos(b/a))

            if duration is None and ecc == 0:
                duration = per/np.pi * np.arcsin(np.sqrt((1 + rp)**2 - b**2) /
                                                 np.sin(np.radians(inc)) / a)

        self.duration = duration
        self.inc = inc
        self.b = b


    @classmethod
    def from_hat11(cls):
        """
        Return a planet with the parameters of HAT-P-11 b.
        """
        return cls(**hat11_params)

    @classmethod
    def from_hat7(cls):
        """
        Return a planet with the parameters of HAT-P-7 b.
        """
        return cls(**hat7_params)

    @classmethod
    def non_transiting(cls):
        """
        Return a planet that does not transit.
        """
        return cls(**null_planet_params)