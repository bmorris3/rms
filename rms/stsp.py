# Licensed under the MIT License - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from glob import glob
import datetime
from threading import Lock
from warnings import warn
import os, subprocess, shutil

import numpy as np
from astropy.io import ascii
from astropy.io.ascii import InconsistentTableError

from .lightcurve import LightCurve
from .exceptions import (OverlappingSpotsWarning, STSPMemoryWarning,
                         STSPFailureWarning)


__all__ = ['STSP', 'clean_up_rms_stsp_dirs']

lock = Lock()

stsp_executable = os.getenv('STSP_EXECUTABLE')

infile_template_l = """#PLANET PROPERTIES
1							; Number of planets -- (if there are more than 1 planet, then the set of 8 planet properties are repeated)
{t0:2.10f}					; T0, epoch         (middle of first transit) in days.
{period:2.10f}				; Planet Period      (days)
{depth:2.10f}				; (Rp/Rs)^2         (Rplanet / Rstar )^ 2
{duration:2.10f}			; Duration (days)   (physical duration of transit, not used)
{b:2.10f}					; Impact parameter  (0= planet cross over equator)
{inclination:2.10f}			; Inclination angle of orbit (90 deg = planet crosses over equator)
{lam:2.10f}					; Lambda of orbit (0 deg = orbital axis along z-axis)
{ecosw:2.10f}			; ecosw
{esinw:2.10f}			; esinw
#STAR PROPERTIES
{rho_s:2.10f} 			; Mean Stellar density (Msun/Rsun^3)
{per_rot:2.10f}			; Stellar Rotation period (days)
4780					; Stellar Temperature
0.31					; Stellar metallicity
{tilt_from_z:2.10f}						; Tilt of the rotation axis of the star down from z-axis (degrees)
{nonlinear_ld}			; Limb darkening (4 coefficients)
{n_ld_rings:d}			; number of rings for limb darkening appoximation
#SPOT PROPERTIES
{n_spots}						; number of spots
{spot_contrast}					; fractional lightness of spots (0.0=total dark, 1.0=same as star)
#LIGHT CURVE
{model_path}			; lightcurve input data file
{start_time:2.10f}		; start time to start fitting the light curve
{lc_duration:2.10f}		; duration of light curve to fit (days)
{real_max:2.10f}		; real maximum of light curve data (corrected for noise), 0 -> use downfrommax
{normalize_oot:d}    	; is light curve flattened (to zero) outside of transits?
#ACTION
l						; l= generate light curve from parameters
{spot_params}
1.00
"""

spot_params_template = """{spot_radius:2.10f}		; spot radius
{spot_theta:2.10f}		; theta
{spot_phi:2.10f}		; phi
"""


def quadratic_to_nonlinear(u1, u2):
    a1 = a3 = 0
    a2 = u1 + 2*u2
    a4 = -u2
    return a1, a2, a3, a4


def _spot_obj_to_params(spot, quiet=False):

    if hasattr(spot, '__len__'):
        non_overlapping_spot_inds = find_overlapping_spots(spot)
        return np.concatenate([[s.r, s.theta, s.phi]
                               for i, s in enumerate(spot)
                               if i in non_overlapping_spot_inds])
    else:
        return np.array([spot.r, spot.theta, spot.phi])


def find_overlapping_spots(spot_list):
    from overlap import find_overlapping_spots as find
    return find(np.array([spot.theta for spot in spot_list]),
                np.array([spot.phi for spot in spot_list]),
                np.array([spot.r for spot in spot_list]))

# def find_overlapping_spots(spot_list, tolerance=1.01, quiet=False):
#     """
#     Find overlapping spots in a list of spot objects.
#
#     Parameters
#     ----------
#     spot_list : list
#     tolerance : float
#     """
#     overlapping_pairs = []
#     spots_with_overlap = []
#     for i in range(len(spot_list)):
#         for j in range(len(spot_list)):
#             if i < j:
#                 sep = np.arccos(np.cos(spot_list[i].theta) *
#                                 np.cos(spot_list[j].theta) +
#                                 np.sin(spot_list[i].theta) *
#                                 np.sin(spot_list[j].theta) *
#                                 np.cos(spot_list[i].phi - spot_list[j].phi))
#                 if sep < tolerance * (spot_list[i].r + spot_list[j].r):
#                     overlapping_pairs.append((i, j))
#
#                     if i not in spots_with_overlap:
#                         spots_with_overlap.append(i)
#                     if j not in spots_with_overlap:
#                         spots_with_overlap.append(j)
#
#     spots_without_overlap = [spot for i, spot in enumerate(spot_list)
#                              if i not in spots_with_overlap]
#     save_these_spot_indices = [i[0] for i in overlapping_pairs]
#     save_these_spots = [spot for i, spot in enumerate(spot_list)
#                         if i in save_these_spot_indices]
#
#     if len(spots_with_overlap) > 0 and not quiet:
#         toss_these_spot_indices = [i[1] for i in overlapping_pairs]
#         toss_these_spots = [spot for i, spot in enumerate(spot_list)
#                             if i in toss_these_spot_indices]
#         warning_message = ('Some spots were overlapping. Tossing one of the two'
#                            ' overlapping spots. \n\nSpots tossed:\n\n' +
#                            '\n'.join(map(str, toss_these_spots)))
#         warn(warning_message, OverlappingSpotsWarning)
#
#     return spots_without_overlap #+ save_these_spots


def get_rms_dirs():
    """
    Get a list of the STSP directories generated by the rms package.
    """
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    rms_stsp_dirs = glob(os.path.join(pkg_dir, '.rms*'))
    return rms_stsp_dirs


def check_for_undeleted_stsp_dirs(n_dirs_threshold=5):
    """
    Check for necessary cleanup warnings.

    The current API writes STSP input files and output files within the rms
    package directory. This function will raise a warning if the number of
    STSP run directories is piling up, so you don't accidentally store lots of
    hidden directories within the rms package directory.
    """
    rms_stsp_dirs = get_rms_dirs()

    if len(rms_stsp_dirs) > n_dirs_threshold:
        warn_message = ("You currently have more than {0} STSP input/output "
                        "directories stored by rms. Consider turning off "
                        "the `keep_dirs` option, and using the "
                        "`rms.clean_up_rms_stsp_dirs` command to delete them "
                        "(and to prevent saving tons of hidden folders within "
                        "the rms package). "
                        "\n\nDirectories to delete:\n\n{1}"
                        .format(n_dirs_threshold, '\n'.join(rms_stsp_dirs)))
        warn(warn_message, STSPMemoryWarning)


def clean_up_rms_stsp_dirs():
    """
    Delete lingering STSP input/output directories created by rms.

    Whenever the STSP object is constructed and `outdir` is not specified, a
    new directory is created within the rms package directory. If you use the
    `keep_dirs` option or an unexpected error occurs while running STSP, those
    hidden directories created by rms can accrue. And nobody wants that. So you
    can use this command to make them go away.
    """
    rms_stsp_dirs = get_rms_dirs()
    for directory in rms_stsp_dirs:
        shutil.rmtree(directory)


class STSP(object):
    """
    Context manager for working with STSP
    """
    def __init__(self, times, star, spot, outdir=None, keep_dir=False,
                 quiet=False):
        """
        Parameters
        ----------
        times : `~astropy.time.Time`
            Time array object
        star : `~rms.Star`
            Parameters for star
        spot : `~rms.Spot` or list of `~rms.Spot` objects
            Spot parameter object(s)
        outdir : str
            Directory to write temporary outputs into
        skip_overlap_check : bool
            If True, skip check that no spots are overlapping.
        """
        self.times = times
        self.star = star
        self.quiet = quiet
        self.spot_params = _spot_obj_to_params(spot, quiet=quiet)
        self.spot_contrast = self.star.spot_contrast

        current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        random_integer = np.random.randint(0, 1e6)

        if outdir is None:
            self.outdir = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          '.rms_tmp_{0}_{1}'
                                          .format(current_time, random_integer)))
        else:
            self.outdir = outdir

        os.makedirs(self.outdir)

        self.model_path = os.path.join(self.outdir, 'model_lc.dat')
        self.keep_dir = keep_dir

    def __enter__(self):
        check_for_undeleted_stsp_dirs()
        return self

    def __exit__(self, *args):
        if not self.keep_dir:
            shutil.rmtree(self.outdir)

    def generate_lightcurve(self, n_ld_rings=40, stsp_exec=None,
                            normalize_oot=False):
        """
        Generate a light curve with STSP.

        Parameters
        ----------
        n_ld_rings : int
            Number of concentric rings to use in the limb-darkening
            approximation
        stsp_exec : str (optional)
            Optionally pass in a path to a different STSP executable with this
            argument.
        normalize_oot : bool
            Normalize the out-of-transit portions of the light curve? Default
            is `False`. Set to `True` when studying spot occultations during
            transits.

        Return
        ------
        lc : `~rms.LightCurve`
            Light curve object with the model from STSP.
        """

        if stsp_exec is None:
            stsp_exec = stsp_executable

        # Normalize light curve to unity
        real_max = 1

        times = self.times.jd
        fluxes = np.ones_like(times)

        np.savetxt(self.model_path, np.vstack([times, fluxes, fluxes]).T,
                   fmt=str('%1.10f'), delimiter='\t', header='stspinputs')

        # Calculate parameters for STSP:
        eccentricity, omega = self.star.planet.ecc, self.star.planet.w
        ecosw = eccentricity*np.cos(np.radians(omega))
        esinw = eccentricity*np.sin(np.radians(omega))
        start_time = times[0]
        lc_duration = times[-1] - times[0] + 1e-6
        nonlinear_ld = quadratic_to_nonlinear(*self.star.u)
        nonlinear_ld_string = ' '.join(map("{0:.5f}".format, nonlinear_ld))

        # get spot parameters sorted out
        spot_params_str = _spot_params_to_string(self.spot_params)

        # Stick those values into the template file

        params_dict = dict(period=self.star.planet.per, ecosw=ecosw,
                esinw=esinw, lam=self.star.planet.lam,
                           tilt_from_z=90-self.star.inc_stellar,
                           start_time=start_time, lc_duration=lc_duration,
                           real_max=real_max, per_rot=self.star.per_rot,
                           rho_s=self.star.rho_s, depth=self.star.planet.rp ** 2,
                           duration=self.star.planet.duration,
                           t0=self.star.planet.t0, b=self.star.planet.b,
                           inclination=self.star.planet.inc,
                           nonlinear_ld=nonlinear_ld_string,
                           n_ld_rings=n_ld_rings,
                           spot_params=spot_params_str[:-1],
                           n_spots=int(len(self.spot_params)/3),
                           model_path=os.path.basename(self.model_path),
                           spot_contrast=self.spot_contrast,
                           normalize_oot=int(normalize_oot))

        in_file_text = infile_template_l.format(**params_dict)

        # Write out the `.in` file
        with open(os.path.join(self.outdir, 'test.in'), 'w') as in_file:
            in_file.write(in_file_text)

        try:
            stdout = subprocess.check_call([stsp_exec, 'test.in'],
                                           cwd=self.outdir)
        except subprocess.CalledProcessError as err:
            if not self.quiet:
                warning_message = ("STSP failed - this could be a result of "
                                   "bad inputs.")
                warn(warning_message, STSPFailureWarning)

        # Read the outputs
        lcout_path = os.path.join(self.outdir, 'test_lcout.txt')
        if not os.path.exists(lcout_path) or os.stat(lcout_path).st_size == 0:
            stsp_times = self.times.jd
            stsp_fluxes = np.ones(len(self.times))
            stsp_flag = 0 * np.ones(len(self.times))

        else:
            try:
                tbl = ascii.read(lcout_path,
                                 format='fast_no_header')
                stsp_times, stsp_fluxes, stsp_flag = tbl['col1'], tbl['col4'].data, tbl['col5']
            except InconsistentTableError:
                stsp_times = self.times.jd
                stsp_fluxes = np.ones(len(self.times)) * np.nan
                stsp_flag = 0 * np.ones(len(self.times))
        return LightCurve(times=stsp_times, fluxes=stsp_fluxes, quarters=stsp_flag)


def _spot_params_to_string(spot_params):
    spot_params_str = ""
    for param_set in np.split(spot_params, len(spot_params)/3):
        spot_params_str += spot_params_template.format(spot_radius=param_set[0],
                                                       spot_theta=param_set[1],
                                                       spot_phi=param_set[2])
    return spot_params_str


