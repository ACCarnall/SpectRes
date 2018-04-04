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

def spectres(new_spec_wavs, old_spec_wavs, spec_fluxes, spec_errs=None):

    """ 
    Function for resampling spectra (and optionally associated uncertainties) onto a new wavelength basis.

    Parameters
    ----------

    new_spec_wavs : numpy.ndarray
        Array containing the new wavelength sampling desired for the spectrum or spectra.

    old_spec_wavs : numpy.ndarray
        1D array containing the current wavelength sampling of the spectrum or spectra.

    spec_fluxes : numpy.ndarray
        Array containing spectral fluxes at the wavelengths specified in old_spec_wavs, last dimension must correspond to the shape of old_spec_wavs.
        Extra dimensions before this may be used to include multiple spectra.

    spec_errs : numpy.ndarray (optional)
        Array of the same shape as spec_fluxes containing uncertainties associated with each spectral flux value.
    
    Returns
    -------

    resampled_fluxes : numpy.ndarray
        Array of resampled flux values, first dimension is the same length as new_spec_wavs, other dimensions are the same as spec_fluxes

    resampled_errs : numpy.ndarray
        Array of uncertainties associated with fluxes in resampled_fluxes. Only returned if spec_errs was specified.

    """

    # Generate arrays of left hand side positions and widths for the old and new bins
    spec_lhs = np.zeros(old_spec_wavs.shape[0])
    spec_widths = np.zeros(old_spec_wavs.shape[0])
    spec_lhs = np.zeros(old_spec_wavs.shape[0])
    spec_lhs[0] = old_spec_wavs[0] - (old_spec_wavs[1] - old_spec_wavs[0])/2
    spec_widths[-1] = (old_spec_wavs[-1] - old_spec_wavs[-2])
    spec_lhs[1:] = (old_spec_wavs[1:] + old_spec_wavs[:-1])/2
    spec_widths[:-1] = spec_lhs[1:] - spec_lhs[:-1]

    filter_lhs = np.zeros(new_spec_wavs.shape[0]+1)
    filter_widths = np.zeros(new_spec_wavs.shape[0])
    filter_lhs[0] = new_spec_wavs[0] - (new_spec_wavs[1] - new_spec_wavs[0])/2
    filter_widths[-1] = (new_spec_wavs[-1] - new_spec_wavs[-2])
    filter_lhs[-1] = new_spec_wavs[-1] + (new_spec_wavs[-1] - new_spec_wavs[-2])/2
    filter_lhs[1:-1] = (new_spec_wavs[1:] + new_spec_wavs[:-1])/2
    filter_widths[:-1] = filter_lhs[1:-1] - filter_lhs[:-2]

    # Check that the range of wavelengths to be resampled_fluxes onto falls within the initial sampling region
    if filter_lhs[0] < spec_lhs[0] or filter_lhs[-1] > spec_lhs[-1]:
        raise ValueError("spectres: The new wavelengths specified must fall within the range of the old wavelength values.")

    #Generate output arrays to be populated
    resampled_fluxes = np.zeros(spec_fluxes[...,0].shape + new_spec_wavs.shape)

    if spec_errs is not None:
        if spec_errs.shape != spec_fluxes.shape:
            raise ValueError("If specified, spec_errs must be the same shape as spec_fluxes.")
        else:
            resampled_fluxes_errs = np.copy(resampled_fluxes)

    start = 0
    stop = 0

    # Calculate the new spectral flux and uncertainty values, loop over the new bins
    for j in range(new_spec_wavs.shape[0]):

        # Find the first old bin which is partially covered by the new bin
        while spec_lhs[start+1] <= filter_lhs[j]:
            start += 1

        # Find the last old bin which is partially covered by the new bin
        while spec_lhs[stop+1] < filter_lhs[j+1]:
            stop += 1

        # If the new bin falls entirely within one old bin the are the same the new flux and new error are the same as for that bin
        if stop == start:

            resampled_fluxes[...,j] = spec_fluxes[...,start]
            if spec_errs is not None:
                resampled_fluxes_errs[...,j] = spec_errs[...,start]

        # Otherwise multiply the first and last old bin widths by P_ij, all the ones in between have P_ij = 1 
        else:

            start_factor = (spec_lhs[start+1] - filter_lhs[j])/(spec_lhs[start+1] - spec_lhs[start])
            end_factor = (filter_lhs[j+1] - spec_lhs[stop])/(spec_lhs[stop+1] - spec_lhs[stop])

            spec_widths[start] *= start_factor
            spec_widths[stop] *= end_factor

            # Populate the resampled_fluxes spectrum and uncertainty arrays
            resampled_fluxes[...,j] = np.sum(spec_widths[start:stop+1]*spec_fluxes[...,start:stop+1], axis=-1)/np.sum(spec_widths[start:stop+1])

            if spec_errs is not None:
                resampled_fluxes_errs[...,j] = np.sqrt(np.sum((spec_widths[start:stop+1]*spec_errs[...,start:stop+1])**2, axis=-1))/np.sum(spec_widths[start:stop+1])
            
            # Put back the old bin widths to their initial values for later use
            spec_widths[start] /= start_factor
            spec_widths[stop] /= end_factor


    # If errors were supplied return the resampled_fluxes spectrum and error arrays
    if spec_errs is not None:
        return resampled_fluxes, resampled_fluxes_errs

    # Otherwise just return the resampled_fluxes spectrum array
    else: 
        return resampled_fluxes


