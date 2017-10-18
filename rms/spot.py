# Licensed under the MIT License - see LICENSE

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.coordinates import Angle

__all__ = ['Spot']


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

    @classmethod
    def at_random_position(cls, radius):
        lat, lon = generate_random_coord()
        return cls(lat, lon, radius)

    @property
    def r(self):
        """Spot radius (alias)"""
        return self.radius

    @property
    def theta(self):
        """Spot ``theta`` [radians] where ``theta = 90 deg - latitude``"""
        return np.pi/2 - self.latitude.to(u.rad).value

    @property
    def phi(self):
        """Spot ``phi`` [radians] where ``phi = longitude``"""
        return self.longitude.to(u.rad).value

    def __repr__(self):
        return "<Spot: lat={0}, lon={1}, rad={2}>".format(self.latitude,
                                                          self.longitude,
                                                          self.radius)

    def plot(self, ax=None, projection='hammer',
             plot_kwargs=dict(marker=',', color='k')):
        """
        Make a simple plot of this spot.

        Parameters
        ----------
        ax : `~matplotlib.pyplot.Axes`
            Axis object
        plot_kwargs : dict
            Keyword arguments to pass to `~matplotlib.pyplot.plot`

        Returns
        -------
        ax : `~matplotlib.pyplot.Axes`
            Updated axis object
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection=projection)

        phi = Angle(self.longitude).wrap_at(np.pi*u.rad)
        theta = (np.pi/2 * u.rad - self.latitude)

        lat, lon = rtp_to_edge(self.radius, phi, theta)
        ax.plot(lon, lat, **plot_kwargs)
        ax.grid()
        return ax


@u.quantity_input(theta=u.rad, phi=u.rad)
def rtp_to_edge(radius, theta, phi, n_points=1000):
    """
    Use the haversine formula to compute the boundary lat/lon coordinates for a
    circular spot centered on ``(theta, phi)``, with radius ``radius` in units
    of the stellar radius.

    Parameters
    ----------
    radius : float
        Spot radius [R_star]
    theta : `~astropy.units.Quantity`
        Spot theta coord (90 deg - latitude)
    phi : `~astropy.units.Quantity`
        Spot phi coord (longitude)
    n_points : int
        Number of points to include in the circle boundary

    Returns
    -------
    lat : `~astropy.units.Quantity`
        Latitudes of spot boundary
    lon : `~astropy.units.Quantity`
        Longitudes of spot boundary
    """
    lat1 = np.pi/2 - theta.to(u.rad).value
    lon1 = phi.to(u.rad).value
    d = radius

    thetas = np.linspace(0, -2*np.pi, n_points)[:, np.newaxis]

    lat = np.arcsin(np.sin(lat1) * np.cos(d) + np.cos(lat1) *
                    np.sin(d) * np.cos(thetas))
    dlon = np.arctan2(np.sin(thetas) * np.sin(d) * np.cos(lat1),
                      np.cos(d) - np.sin(lat1) * np.sin(lat))
    lon = ((lon1 - dlon + np.pi) % (2*np.pi)) - np.pi
    return lat*u.rad, lon*u.rad


def generate_random_coord(n=1):
    """
    Generate random latitude/longitude pairs, of length ``n``.
    """
    random_longitude = 2*np.pi*np.random.rand(n)
    random_z = 2*np.random.rand(n) - 1
    random_latitude = np.arcsin(random_z)
    result = np.vstack([random_latitude, random_longitude]).T * u.rad

    if result.shape[0] == 1:
        return result[0]
    return result
