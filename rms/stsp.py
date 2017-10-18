# Licensed under the MIT License - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from glob import glob
import datetime
import os, subprocess, shutil, time
import numpy as np
from astropy.io import ascii
from .lightcurve import LightCurve

from threading import Lock

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
0						; is light curve flattened (to zero) outside of transits?
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


def _clean_up(require_input=False):
    paths_to_clean = glob(os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                       '.rms_tmp_*')))
    if require_input:
        user_input = input("Delete following paths [y]/n: \n" +
                           '\n'.join(paths_to_clean))
        if not user_input.lower() == 'n':
            for directory in paths_to_clean:
                shutil.rmtree(directory)
    else:
        for directory in paths_to_clean:
            shutil.rmtree(directory)


def _spot_obj_to_params(spot):
    if hasattr(spot, '__len__'):
        return np.concatenate([[s.r, s.theta, s.phi] for s in spot])
    else:
        return np.array([spot.r, spot.theta, spot.phi])


class STSP(object):
    """
    Context manager for working with STSP
    """
    def __init__(self, times, star, spot, outdir=None, keep_dir=False):
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
        """
        self.times = times
        self.star = star
        self.spot_params = _spot_obj_to_params(spot)
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
        return self

    def __exit__(self, *args):
        if not self.keep_dir:
            shutil.rmtree(self.outdir)

    def safe_clean_up(self):
        paths_to_delete = ['model_lc.dat', 'test.in', 'xyzdetail.txt',
                           'test_lcout.txt', 'test_errstsp.txt']
        for path in paths_to_delete:
            abspath = os.path.join(self.outdir, path)
            if os.path.exists(abspath):
                os.remove(abspath)

    def generate_lightcurve(self, n_ld_rings=40, stsp_exec=None):
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

        Return
        ------
        lc : `~rms.LightCurve`
            Light curve object with the model from STSP.
        """

        if stsp_exec is None:
            stsp_exec = stsp_executable

        # Normalize light curve to unity
        real_max = 1

        n_transits = np.rint(np.median((self.star.t0 -
                                        self.times.jd) /
                                       self.star.per))

        times = self.times.jd
        fluxes = np.ones_like(times)

        np.savetxt(self.model_path,
                   np.vstack([times, fluxes,
                              fluxes]).T,
                   fmt=str('%1.10f'), delimiter='\t', header='stspinputs')

        # Calculate parameters for STSP:
        eccentricity, omega = self.star.ecc, self.star.w
        ecosw = eccentricity*np.cos(np.radians(omega))
        esinw = eccentricity*np.sin(np.radians(omega))
        start_time = times[0]
        lc_duration = times[-1] - times[0]
        nonlinear_ld = quadratic_to_nonlinear(*self.star.u)
        nonlinear_ld_string = ' '.join(map("{0:.5f}".format, nonlinear_ld))

        # get spot parameters sorted out
        spot_params_str = _spot_params_to_string(self.spot_params)

        # Stick those values into the template file

        params_dict = dict(period=self.star.per, ecosw=ecosw,
                           esinw=esinw, lam=self.star.lam,
                           tilt_from_z=90-self.star.inc_stellar,
                           start_time=start_time, lc_duration=lc_duration,
                           real_max=real_max, per_rot=self.star.per_rot,
                           rho_s=1.0, depth=self.star.rp ** 2,
                           duration=self.star.duration,
                           t0=self.star.t0, b=self.star.b,
                           inclination=self.star.inc,
                           nonlinear_ld=nonlinear_ld_string,
                           n_ld_rings=n_ld_rings,
                           spot_params=spot_params_str[:-1],
                           n_spots=int(len(self.spot_params)/3),
                           model_path=os.path.basename(self.model_path),
                           spot_contrast=self.spot_contrast)

        in_file_text = infile_template_l.format(**params_dict)

        # Write out the `.in` file
        with open(os.path.join(self.outdir, 'test.in'), 'w') as in_file:
            in_file.write(in_file_text)

        try:
            stdout = subprocess.check_output([stsp_exec, 'test.in'], cwd=self.outdir)
        except subprocess.CalledProcessError as err:
            pass  # print("Failed. Error:", err.output, err.stderr, err.stdout)

        time.sleep(0.01)
        # Read the outputs
        if os.stat(os.path.join(self.outdir, 'test_lcout.txt')).st_size == 0:
            stsp_times = self.times.jd
            stsp_fluxes = np.ones_like(self.times)
            stsp_flag = 0 * np.ones_like(self.times)

        else:
            tbl = ascii.read(os.path.join(self.outdir, 'test_lcout.txt'),
                             format='fast_no_header')
            stsp_times, stsp_fluxes, stsp_flag = tbl['col1'], tbl['col4'], tbl['col5']
        return LightCurve(times=stsp_times, fluxes=stsp_fluxes.data, quarters=stsp_flag)


def _spot_params_to_string(spot_params):
    spot_params_str = ""
    for param_set in np.split(spot_params, len(spot_params)/3):
        spot_params_str += spot_params_template.format(spot_radius=param_set[0],
                                                       spot_theta=param_set[1],
                                                       spot_phi=param_set[2])
    return spot_params_str


