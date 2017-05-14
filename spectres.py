"""
SpectRes: A fast spectral resampling function.
Copyright (C) 2017  A. C. Carnall

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


import numpy as np 
import sys


# Function for calculating the left hand side (lhs) positions and widths of the spectral bins from their central wavelengths.
def make_bins(wavelengths, make_rhs="False"):
    bin_widths = np.zeros(wavelengths.shape[0])

    # This option makes the final entry in the left hand sides array the right hand side of the final bin
    if make_rhs == "True":
        bin_lhs = np.zeros(wavelengths.shape[0]+1)
        #The first lhs position is assumed to be as far from the first central wavelength as the rhs of the first bin.
        bin_lhs[0] = wavelengths[0] - (wavelengths[1]-wavelengths[0])/2
        bin_widths[-1] = (wavelengths[-1] - wavelengths[-2])
        bin_lhs[-1] = wavelengths[-1] + (wavelengths[-1]-wavelengths[-2])/2
        bin_lhs[1:-1] = (wavelengths[1:] + wavelengths[:-1])/2
        bin_widths[:-1] = bin_lhs[1:-1]-bin_lhs[:-2]

    # Otherwise just return the lhs positions of each bin
    else:
        bin_lhs = np.zeros(wavelengths.shape[0])
        bin_lhs[0] = wavelengths[0] - (wavelengths[1]-wavelengths[0])/2
        bin_widths[-1] = (wavelengths[-1] - wavelengths[-2])
        bin_lhs[1:] = (wavelengths[1:] + wavelengths[:-1])/2
        bin_widths[:-1] = bin_lhs[1:]-bin_lhs[:-1]

    return bin_lhs, bin_widths



# Function for performing spectral resampling on a spectrum or array of spectra.
def spectres(spec_wavs, spec_fluxes, resampling, spec_errs=None):

    # Generate arrays of left hand side positions and widths for the old and new bins
    filter_lhs, filter_widths = make_bins(resampling, make_rhs="True")
    spec_lhs, spec_widths = make_bins(spec_wavs)


    # Check that the range of wavelengths to be resampled onto falls within the initial sampling region
    if filter_lhs[0] < spec_lhs[0] or filter_lhs[-1] > spec_lhs[-1]:
        print "Spec_lhs, filter_lhs, filter_rhs, spec_rhs ", spec_lhs[0], filter_lhs[0], filter_lhs[-1], spec_lhs[-1]
        sys.exit("spectres was passed a spectrum which did not cover the full wavelength range of the specified filter curve.")
    

    #Generate output arrays to be populated
    if spec_fluxes.ndim == 1:
        resampled = np.zeros((resampling.shape[0]))

    elif spec_fluxes.ndim == 2:
        resampled = np.zeros((len(resampling), spec_fluxes.shape[1]))

    if spec_errs is not None:
        if spec_errs.shape != spec_fluxes.shape:
            sys.exit("If specified, spec_errs must be the same shape as spec_fluxes.")
        else:
            resampled_errs = np.copy(resampled)

    start = 0
    stop = 0

    # Calculate the new spectral flux and uncertainty values, loop over the new bins
    for j in range(len(filter_lhs)-1):

        # Find the first old bin which is partially covered by the new bin
        while spec_lhs[start+1] <= filter_lhs[j]:
            start += 1

        # Find the last old bin which is partially covered by the new bin
        while spec_lhs[stop+1] < filter_lhs[j+1]:
            stop += 1

        if spec_fluxes.ndim == 1:

            # If the new bin falls entirely within one old bin the are the same the new flux and new error are the same as for that bin
            if stop == start:

                resampled[j] = spec_fluxes[start]
                if spec_errs is not None:
                    resampled_errs[j] = spec_errs[start]

            # Otherwise multiply the first and last old bin widths by P_ij, all the ones in between have P_ij = 1 
            else:

                start_factor = (spec_lhs[start+1] - filter_lhs[j])/(spec_lhs[start+1] - spec_lhs[start])
                end_factor = (filter_lhs[j+1] - spec_lhs[stop])/(spec_lhs[stop+1] - spec_lhs[stop])

                spec_widths[start] *= start_factor
                spec_widths[stop] *= end_factor

                # Populate the resampled spectrum and uncertainty arrays
                resampled[j] = np.sum(spec_widths[start:stop+1]*spec_fluxes[start:stop+1])/np.sum(spec_widths[start:stop+1])

                if spec_errs is not None:
                    resampled_errs[j] = np.sqrt(np.sum((spec_widths[start:stop+1]*spec_errs[start:stop+1])**2))/np.sum(spec_widths[start:stop+1])
                
                # Put back the old bin widths to their initial values for later use
                spec_widths[start] /= start_factor
                spec_widths[stop] /= end_factor


        # The same as above, except operates on each row of the array, resampling all of the input models
        elif spec_fluxes.ndim == 2:

            if stop == start:

                resampled[j, :] = spec_fluxes[start, :]
                if spec_errs is not None:
                    resampled_errs[j, :] = spec_errs[start, :]

            else:

                start_factor = (spec_lhs[start+1] - filter_lhs[j])/(spec_lhs[start+1] - spec_lhs[start])
                end_factor = (filter_lhs[j+1] - spec_lhs[stop])/(spec_lhs[stop+1] - spec_lhs[stop])

                spec_widths[start] *= start_factor
                spec_widths[stop] *= end_factor

                resampled[j, :] = np.sum(np.expand_dims(spec_widths[start:stop+1], axis=1)*spec_fluxes[start:stop+1, :], axis=0)/np.sum(spec_widths[start:stop+1])
                
                if spec_errs is not None:
                    resampled_errs[j, :] = np.sqrt(np.sum((np.expand_dims(spec_widths[start:stop+1], axis=1)*spec_errs[start:stop+1])**2, axis=0))/np.sum(spec_widths[start:stop+1])
                
                spec_widths[start] /= start_factor
                spec_widths[stop] /= end_factor


    # If errors were supplied return the resampled spectrum and error arrays
    if spec_errs is not None:
        return resampled, resampled_errs

    # Otherwise just return the resampled spectrum array
    else: 
        return resampled


