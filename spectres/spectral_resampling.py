from __future__ import print_function, division, absolute_import
import warnings

import numpy as np


def make_bins(wavs):
    """ Given a series of wavelength points, find the edges and widths
    of corresponding wavelength bins. """
    edges = np.zeros(wavs.shape[0]+1)
    widths = np.zeros(wavs.shape[0])
    edges[0] = wavs[0] - (wavs[1] - wavs[0])/2
    widths[-1] = (wavs[-1] - wavs[-2])
    edges[-1] = wavs[-1] + (wavs[-1] - wavs[-2])/2
    edges[1:-1] = (wavs[1:] + wavs[:-1])/2
    widths[:-1] = edges[1:-1] - edges[:-2]

    return edges, widths


def spectres(new_wavs, spec_wavs, spec_fluxes, spec_errs=None, fill=None,
             verbose=True):

    """
    Function for resampling spectra (and optionally associated
    uncertainties) onto a new wavelength basis.

    Parameters
    ----------

    new_wavs : numpy.ndarray
        Array containing the new wavelength sampling desired for the
        spectrum or spectra.

    spec_wavs : numpy.ndarray
        1D array containing the current wavelength sampling of the
        spectrum or spectra.

    spec_fluxes : numpy.ndarray
        Array containing spectral fluxes at the wavelengths specified in
        spec_wavs, last dimension must correspond to the shape of
        spec_wavs. Extra dimensions before this may be used to include
        multiple spectra.

    spec_errs : numpy.ndarray (optional)
        Array of the same shape as spec_fluxes containing uncertainties
        associated with each spectral flux value.

    fill : float (optional)
        Where new_wavs extends outside the wavelength range in spec_wavs
        this value will be used as a filler in new_fluxes and new_errs.

    verbose : bool (optional)
        Setting verbose to False will suppress the default warning about
        new_wavs extending outside spec_wavs and "fill" being used.

    Returns
    -------

    new_fluxes : numpy.ndarray
        Array of resampled flux values, last dimension is the same
        length as new_wavs, other dimensions are the same as
        spec_fluxes.

    new_errs : numpy.ndarray
        Array of uncertainties associated with fluxes in new_fluxes.
        Only returned if spec_errs was specified.
    """

    # Rename the input variables for clarity within the function.
    old_wavs = spec_wavs
    old_fluxes = spec_fluxes
    old_errs = spec_errs

    # Make arrays of edge positions and widths for the old and new bins

    old_edges, old_widths = make_bins(old_wavs)
    new_edges, new_widths = make_bins(new_wavs)

    if old_errs is not None and old_errs.shape != old_fluxes.shape:
        raise ValueError("If specified, spec_errs must be the same shape "
                         "as spec_fluxes.")

    # Find first and last old bins partially covered by each new bin
    last = old_wavs.shape[0] - 1
    start = np.clip(np.searchsorted(old_edges, new_edges[:-1], side="right") - 1,
                    0, last)
    stop = np.clip(np.searchsorted(old_edges, new_edges[1:], side="left") - 1,
                   0, last)

    # Widths of the parts of the first and last old bins which are covered
    start_width = old_edges[start+1] - new_edges[:-1]
    end_width = new_edges[1:] - old_edges[stop]

    # Indices for summing over old bins fully covered by each new bin
    has_middle = (stop - start > 1)
    sum_inds = np.minimum(np.c_[start+1, stop].ravel(), last)

    # If new bin is fully inside an old bin start and stop are equal
    single = (stop == start)

    # Bins where new_wavs extends outside of spec_wavs get filler values
    outside = (new_edges[:-1] < old_edges[0]) | (new_edges[1:] > old_edges[-1])

    if verbose and (outside[0] or outside[-1]):
        warnings.warn(
            "Spectres: new_wavs contains values outside the range "
            "in spec_wavs, new_fluxes and new_errs will be filled "
            "with the value set in the 'fill' keyword argument "
            "(by default 0).",
            category=RuntimeWarning,
        )

    def resample(vals, power):
        """ Sum (width*vals)**power over each new bin, then normalise. """
        middle = np.add.reduceat((old_widths*vals)**power, sum_inds,
                                 axis=-1)[..., ::2]*has_middle

        total = (middle + (start_width*vals[..., start])**power
                 + (end_width*vals[..., stop])**power)**(1./power)

        new_vals = total/new_widths
        new_vals[..., single] = vals[..., start[single]]
        new_vals[..., outside] = fill
        return new_vals

    # Populate new_fluxes spectrum and uncertainty arrays
    new_fluxes = resample(old_fluxes, 1)

    if old_errs is not None:
        new_errs = resample(old_errs, 2)

    # If errors were supplied return both new_fluxes and new_errs.
    if old_errs is not None:
        return new_fluxes, new_errs

    # Otherwise just return the new_fluxes spectrum array
    else:
        return new_fluxes
