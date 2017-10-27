# Licensed under the MIT License - see LICENSE.rst

from astropy.utils.exceptions import AstropyUserWarning

__all__ = ['OverlappingSpotsWarning', 'STSPMemoryWarning', 'STSPFailureWarning']


class OverlappingSpotsWarning(AstropyUserWarning):
    """
    Warning for when a user submits spots that would be overlapping.
    """
    pass


class STSPMemoryWarning(AstropyUserWarning):
    """
    Warning for when a user is accruing lots of STSP input/output directories.
    """
    pass

class STSPFailureWarning(AstropyUserWarning):
    """
    Warning for when STSP seg-faults or returns an initialization error
    """
    pass
