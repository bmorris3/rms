# Licensed under the MIT License - see LICENSE.rst

from astropy.utils.exceptions import AstropyUserWarning

__all__ = ['OverlappingSpotsWarning']

class OverlappingSpotsWarning(AstropyUserWarning):
    """
    Warning for when a user submits spots that would be overlapping.
    """
    pass