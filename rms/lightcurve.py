# Licensed under the MIT License - see LICENSE.rst
"""
Methods for taking the raw light curves from MAST and producing cleaned light
curves.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import os
import numpy as np
import matplotlib.pyplot as plt
import shutil
#import batman


def generate_lc_depth(times, depth, transit_params):
    """
    Generate a model transit light curve.

    Parameters
    ----------
    times : `~numpy.ndarray`
        Times in JD
    depth : float
        Set depth independently from the setting in `transit_params`
    transit_params : `~batman.TransitParams`
        Transit light curve parameters

    Returns
    -------

    """
    exp_time = (1*u.min).to(u.day).value

    transit_params.rp = np.sqrt(depth)

    m = batman.TransitModel(transit_params, times, supersample_factor=7,
                            exp_time=exp_time)
    model_flux = m.light_curve(transit_params)
    return model_flux


class LightCurve(object):
    """
    Container object for light curves.
    """
    def __init__(self, times=None, fluxes=None, errors=None, quarters=None,
                 name=None):
        """
        Parameters
        ----------
        times : `~numpy.ndarray`
            Times in JD
        fluxes : `~numpy.ndarray`
            Fluxes (normalized or not)
        errors : `~numpy.ndarray`
            Uncertainties on the fluxes
        quarters : `~numpy.ndarray` (optional)
            Kepler Quarter for each flux
        name : str
            Name this light curve (optional)
        """
        # if len(times) < 1:
        #    raise ValueError("Input `times` have no length.")

        if isinstance(times[0], Time) and isinstance(times, np.ndarray):
            times = Time(times)
        elif not isinstance(times, Time):
            times = Time(times, format='jd')

        self.times = times
        self.fluxes = fluxes
        if self.times is not None and errors is None:
            errors = np.zeros_like(self.fluxes) - 1
        self.errors = errors
        if self.times is not None and quarters is None:
            quarters = np.zeros_like(self.fluxes) - 1
        self.quarters = quarters
        self.name = name

    def phases(self, params):
        phase = ((self.times.jd - params.t0) % params.per)/params.per
        phase[phase > 0.5] -= 1.0
        return phase

    def plot(self, star=None, ax=None, show=True,
             phase=None, plot_kwargs={'color':'k', 'lw':0, 'marker':'.'}):
        """
        Plot light curve.

        Parameters
        ----------
        star : `~rms.Star` (optional)
            Star parameters. Required if `phase` is `True`.
        ax : `~matplotlib.axes.Axes` (optional)
            Axis to make plot on top of
        show : bool
            If `True`, call `matplotlib.pyplot.show` after plot is made
        phase : bool
            If `True`, map times in JD to orbital phases, which requires
            that `transit_params` be input also.
        plot_kwargs : dict
            Keyword arguments to pass to `~matplotlib` calls.
        """
        if ax is None:
            ax = plt.gca()

        if star is not None and phase is None:
            phase = True

        if phase:
            x = (self.times.jd - star.t0)/star.per_rot % 1
            first_half = x < 0.5
            second_half = x >= 0.5
            ax.plot(x[second_half] - 1, self.fluxes[second_half], '.', color='gray', alpha=0.5)
            ax.plot(x[first_half] + 1, self.fluxes[first_half], '.', color='gray', alpha=0.5)
            ax.plot(x, self.fluxes, **plot_kwargs)

        else:
            ax.plot(self.times.jd, self.fluxes, **plot_kwargs)

        ax.set(xlabel='Time' if not phase else 'Phase', ylabel='Flux')

        if self.name is not None:
            ax.set_title(self.name)
        if show:
            plt.show()
        return ax

    def save_to(self, path, overwrite=False, for_stsp=False):
        """
        Save times, fluxes, errors to new directory ``dirname`` in ``path``
        """
        dirname = self.name
        output_path = os.path.join(path, dirname)
        self.times = Time(self.times)

        if not for_stsp:
            if os.path.exists(output_path) and overwrite:
                shutil.rmtree(output_path)

            if not os.path.exists(output_path):
                os.mkdir(output_path)
                for attr in ['times_jd', 'fluxes', 'errors', 'quarters']:
                    np.savetxt(os.path.join(path, dirname,
                                            '{0}.txt'.format(attr)),
                               getattr(self, attr))

        else:
            if not os.path.exists(output_path) or overwrite:
                attrs = ['times_jd', 'fluxes', 'errors']
                output_array = np.zeros((len(self.fluxes), len(attrs)),
                                        dtype=float)
                for i, attr in enumerate(attrs):
                    output_array[:, i] = getattr(self, attr)
                np.savetxt(os.path.join(path, dirname+'.txt'), output_array)

    @classmethod
    def from_raw_fits(cls, fits_paths, name=None):
        """
        Load FITS files downloaded from MAST into the `LightCurve` object.

        Parameters
        ----------
        fits_paths : list
            List of paths to FITS files to read in
        name : str (optional)
            Name of light curve

        Returns
        -------
        lc : `LightCurve`
            The light curve for the data in the fits files.
        """
        fluxes = []
        errors = []
        times = []
        quarter = []

        # Manual on times: http://archive.stsci.edu/kepler/manuals/archive_manual.htm

        for path in fits_paths:
            data = fits.getdata(path)
            header = fits.getheader(path)
            timslice = fits.open(path)[1].header['TIMSLICE']
            time_slice_correction = (0.25 + 0.62*(5.0 - timslice))/86400
            times.append(data['TIME'] + 2454833.0)# - data['TIMECORR'] + time_slice_correction)
            errors.append(data['SAP_FLUX_ERR'])
            fluxes.append(data['SAP_FLUX'])
            quarter.append(len(data['TIME'])*[header['QUARTER']])

        times, fluxes, errors, quarter = [np.concatenate(i)
                                          for i in [times, fluxes,
                                                    errors, quarter]]

        mask_nans = np.zeros_like(fluxes).astype(bool)
        for attr in [times, fluxes, errors]:
            mask_nans |= np.isnan(attr)

        times, fluxes, errors, quarter = [attr[-mask_nans]
                                           for attr in [times, fluxes, errors, quarter]]

        return LightCurve(times, fluxes, errors, quarters=quarter, name=name)

    @classmethod
    def from_dir(cls, path, for_stsp=False):
        """Load light curve from numpy save files in ``dir``"""
        if not for_stsp:
            times, fluxes, errors, quarters = [np.loadtxt(os.path.join(path, '{0}.txt'.format(attr)))
                                               for attr in ['times_jd', 'fluxes', 'errors', 'quarters']]
        else:
            quarters = None
            times, fluxes, errors = np.loadtxt(path, unpack=True)

        if os.sep in path:
            name = path.split(os.sep)[-1]
        else:
            name = path

        if name.endswith('.txt'):
            name = name[:-4]

        return cls(times, fluxes, errors, quarters=quarters, name=name)

    def normalize_each_quarter(self, rename=None, polynomial_order=2,
                               plots=False):
        """
        Use polynomial fit to each quarter to normalize the data.

        Parameters
        ----------
        rename : str (optional)
            New name of the light curve after normalization
        polynomial_order : int (optional)
            Order of polynomial to fit to the out-of-transit fluxes. Default
            is 2.
        plots : bool (optional)
            Show diagnostic plots after normalization.
        """
        quarter_inds = list(set(self.quarters))
        quarter_masks = [quarter == self.quarters for quarter in quarter_inds]

        for quarter_mask in quarter_masks:

            polynomial = np.polyfit(self.times[quarter_mask].jd,
                                    self.fluxes[quarter_mask], polynomial_order)
            scaling_term = np.polyval(polynomial, self.times[quarter_mask].jd)
            self.fluxes[quarter_mask] /= scaling_term
            self.errors[quarter_mask] /= scaling_term

            if plots:
                plt.plot(self.times[quarter_mask], self.fluxes[quarter_mask])
                plt.show()

        if rename is not None:
            self.name = rename

    def delete_outliers(self):

        d = np.diff(self.fluxes)
        spikey = np.abs(d - np.median(d)) > 2.5*np.std(d)
        neighboring_spikes = spikey[1:] & spikey[:-1]
        opposite_signs = np.sign(d[1:]) != np.sign(d[:-1])
        outliers = np.argwhere(neighboring_spikes & opposite_signs) + 1
        #print('number bad fluxes: {0}'.format(len(outliers)))

        self.times = Time(np.delete(self.times.jd, outliers), format='jd')
        self.fluxes = np.delete(self.fluxes, outliers)
        self.errors = np.delete(self.errors, outliers)
        self.quarters = np.delete(self.quarters, outliers)

    def mask_out_of_transit(self, params, oot_duration_fraction=0.25,
                            flip=False):
        """
        Mask out the out-of-transit light curve based on transit parameters

        Parameters
        ----------
        params : `~batman.TransitParams`
            Transit light curve parameters. Requires that `params.duration`
            is defined.
        oot_duration_fraction : float (optional)
            Fluxes from what fraction of a transit duration of the
            out-of-transit light curve should be included in the mask?
        flip : bool (optional)
            If `True`, mask in-transit rather than out-of-transit.

        Returns
        -------
        d : dict
            Inputs for a new `LightCurve` object with the mask applied.
        """
        # Fraction of one duration to capture out of transit

        phased = (self.times.jd - params.t0) % params.per
        near_transit = ((phased < params.duration*(0.5 + oot_duration_fraction)) |
                        (phased > params.per - params.duration*(0.5 + oot_duration_fraction)))
        if flip:
            near_transit = -near_transit
        sort_by_time = np.argsort(self.times[near_transit].jd)
        return dict(times=self.times[near_transit][sort_by_time],
                    fluxes=self.fluxes[near_transit][sort_by_time],
                    errors=self.errors[near_transit][sort_by_time],
                    quarters=self.quarters[near_transit][sort_by_time])

    def mask_in_transit(self, params, oot_duration_fraction=0.25):
        """
        Mask out the in-transit light curve based on transit parameters

        Parameters
        ----------
        params : `~batman.TransitParams`
            Transit light curve parameters. Requires that `params.duration`
            is defined.
        oot_duration_fraction : float (optional)
            Fluxes from what fraction of a transit duration of the
            out-of-transit light curve should be included in the mask?

        Returns
        -------
        d : dict
            Inputs for a new `LightCurve` object with the mask applied.
        """
        return self.mask_out_of_transit(params, flip=True,
                                        oot_duration_fraction=oot_duration_fraction)

    def get_transit_light_curves(self, params, plots=False):
        """
        For a light curve with transits only (i.e. like one returned by
        `LightCurve.mask_out_of_transit`), split up the transits into their
        own light curves, return a list of `TransitLightCurve` objects.

        Parameters
        ----------
        params : `~batman.TransitParams`
            Transit light curve parameters

        plots : bool
            Make diagnostic plots.

        Returns
        -------
        transit_light_curves : list
            List of `TransitLightCurve` objects
        """
        time_diffs = np.diff(sorted(self.times.jd))
        diff_between_transits = params.per/2.
        split_inds = np.argwhere(time_diffs > diff_between_transits) + 1

        if len(split_inds) > 1:

            split_ind_pairs = [[0, split_inds[0][0]]]
            split_ind_pairs.extend([[split_inds[i][0], split_inds[i+1][0]]
                                     for i in range(len(split_inds)-1)])
            split_ind_pairs.extend([[split_inds[-1], len(self.times)]])

            transit_light_curves = []
            counter = -1
            for start_ind, end_ind in split_ind_pairs:
                counter += 1
                if plots:
                    plt.plot(self.times.jd[start_ind:end_ind],
                             self.fluxes[start_ind:end_ind], '.-')

                parameters = dict(times=self.times[start_ind:end_ind],
                                  fluxes=self.fluxes[start_ind:end_ind],
                                  errors=self.errors[start_ind:end_ind],
                                  quarters=self.quarters[start_ind:end_ind],
                                  name=counter)
                transit_light_curves.append(TransitLightCurve(**parameters))
            if plots:
                plt.show()
        else:
            transit_light_curves = []

        return transit_light_curves

    def get_available_quarters(self):
        """
        Get which quarters are available in this `LightCurve`

        Returns
        -------
        qs : list
            List of unique quarters available.
        """
        return list(set(self.quarters))

    def get_quarter(self, quarter):
        """
        Get a copy of the data from within `LightCurve` during one Kepler
        quarter.

        Parameters
        ----------
        quarter : int
            Kepler Quarter

        Returns
        -------
        lc : `LightCurve`
            Light curve from one Kepler Quarter
        """
        this_quarter = self.quarters == quarter
        return LightCurve(times=self.times[this_quarter],
                          fluxes=self.fluxes[this_quarter],
                          errors=self.errors[this_quarter],
                          quarters=self.quarters[this_quarter],
                          name=self.name + '_quarter_{0}'.format(quarter))

    @property
    def times_jd(self):
        """
        Get the times in this light curve in JD.

        Returns
        -------
        t_jd : `~numpy.ndarray`
            Julian dates.
        """
        return self.times.jd

    def split_at_index(self, index):
        """
        Split the light curve into two light curves, at ``index``
        """
        return (LightCurve(times=self.times[:index], fluxes=self.fluxes[:index], 
                           errors=self.errors[:index], quarters=self.quarters[:index], 
                           name=self.name),
                LightCurve(times=self.times[index:], fluxes=self.fluxes[index:], 
                           errors=self.errors[index:], quarters=self.quarters[index:], 
                           name=self.name))

    def transit_model(self, transit_params, short_cadence=True):
        # (1 * u.min).to(u.day).value
        if short_cadence:
            exp_time = (1 * u.min).to(u.day).value #(6.019802903 * 10 * u.s).to(u.day).value
            supersample = 10
        else:
            exp_time = (6.019802903 * 10 * 30 * u.s).to(u.day).value
            supersample = 10

        m = batman.TransitModel(transit_params, self.times.jd,
                                supersample_factor=supersample,
                                exp_time=exp_time)
        model_flux = m.light_curve(transit_params)
        return model_flux
