import os
import sys
from itertools import product
sys.path.insert(0, os.path.abspath('../'))

import h5py
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.utils.console import ProgressBar
from scipy.stats import binned_statistic

from rms import STSP, Star, Spot, Planet

hdf5_file = h5py.File('results.hdf5', 'w')

t0 = 2454954.357463
rotation_period = 25
spot_contrast = 0.2
n_transits = 150
mean_latitudes = np.arange(15, 95, 10)
inc_stellars = np.arange(-60, 70, 10)
lams = [182.5, -132.6, 155]

transit_times = Time(t0, format='jd') + np.arange(-3/24, 3/24, 1/60/24) * u.day

planet = Planet.from_hat7()

with ProgressBar(len(lams) * len(mean_latitudes) * len(inc_stellars)) as bar:
    for lam in lams:
        lam_group = hdf5_file.create_group('lam{0:.2f}'.format(lam))
        lam_group.attrs['lam'] = lam

        for mean_latitude, inc_stellar in product(mean_latitudes, inc_stellars):

            group = lam_group.create_group('latitude{0:.2f}_incstellar{1:.2f}'
                                           .format(mean_latitude, inc_stellar))
            group.attrs['mean_latitude'] = mean_latitude
            group.attrs['inc_stellar'] = inc_stellar

            star = Star(planet=planet, rotation_period=rotation_period,
                        inc_stellar=inc_stellar, spot_contrast=spot_contrast)

            fluxes = np.zeros((len(transit_times), n_transits))
            times = np.zeros((len(transit_times), n_transits))

            no_transit_spot = Spot(radius=0.0001, latitude=-88*u.deg, longitude=0*u.deg)

            i = 0
            while i < n_transits:

                # Draw spots from the sunspot distribution
                n_spots = 50
                spots = [Spot.from_sunspot_distribution(mean_latitude=mean_latitude)
                         for i in range(n_spots)]

                # Draw observing times from uniform distribution
                random_times = transit_times + 1 * (np.random.rand() - 0.5) * u.min

                with STSP(random_times, star, spots, quiet=True) as stsp:
                    lc = stsp.generate_lightcurve(n_ld_rings=100)

                with STSP(random_times, star, no_transit_spot, quiet=True) as stsp:
                    no_spot_lc = stsp.generate_lightcurve(n_ld_rings=100)

                # If STSP ran successfully:
                f_test = lc.fluxes if not hasattr(lc.fluxes, 'value') else lc.fluxes.value

                if not np.any(np.isnan(lc.fluxes)):
                    oot = lc.mask_out_of_transit(star, flip=True)
                    oot_fluxes = oot['fluxes']
                    median_oot_flux = np.median(oot_fluxes)

                    min_length = min([len(lc.fluxes), len(no_spot_lc.fluxes)])
                    if len(lc.fluxes) == len(no_spot_lc.fluxes):
                        fluxes[:, i] = lc.fluxes / median_oot_flux - no_spot_lc.fluxes
                        times[:, i] = random_times.jd

                        i += 1

            group.create_dataset('times', data=times.ravel(), compression="gzip")
            group.create_dataset('fluxes', data=fluxes.ravel(), compression="gzip")

            t, f = times.ravel(), fluxes.ravel()
            f = f[np.argsort(t)]
            t = t[np.argsort(t)]

            n_oneminbins = t.ptp() * 24 * 60
            bs = binned_statistic(t, f, statistic='mean', bins=int(n_oneminbins))
            bincenters = 0.5*(bs.bin_edges[1:] + bs.bin_edges[:-1])
            plt.plot(t, 1e6*f, '.')
            plt.plot(bincenters, 1e6*bs.statistic, 'r', lw=4)
            plt.title(r'$i_s = {0}$, $\ell = {1}^\circ$'.format(inc_stellar, mean_latitude))
            plt.savefig('plots/is{0:.2f}_lat{1:.2f}.png'.format(inc_stellar, mean_latitude),
                        bbox_inches='tight')
            plt.close()

            hdf5_file.flush()
            bar.update()

hdf5_file.close()