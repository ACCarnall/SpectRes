from __future__ import print_function, division, absolute_import
import warnings
import numpy as np


def spectres(new_spec_wavs,
             old_spec_wavs,
             spec_fluxes,
             spec_errs=None,
             fill_value=0.,
             silence=False):
    """
    Function for resampling spectra (and optionally associated
    uncertainties) onto a new wavelength basis.

    Parameters
    ----------
    new_spec_wavs : numpy.ndarray
        Array containing the new wavelength sampling desired for the
        spectrum or spectra.

    old_spec_wavs : numpy.ndarray
        1D array containing the current wavelength sampling of the
        spectrum or spectra.

    spec_fluxes : numpy.ndarray
        Array containing spectral fluxes at the wavelengths specified in
        old_spec_wavs, last dimension must correspond to the shape of
        old_spec_wavs. Extra dimensions before this may be used to
        include multiple spectra.

    spec_errs : numpy.ndarray (optional)
        Array of the same shape as spec_fluxes containing uncertainties
        associated with each spectral flux value.

    fill_value : float (optional, default is 0.)
        Values to be filled for the (new) bins lying outside of the input
        wavelength range.

    silence : tuple (optional, default is False)
        Set to suppress warnings.

    Returns
    -------
    res_fluxes : numpy.ndarray
        Array of resampled flux values, first dimension is the same
        length as new_spec_wavs, other dimensions are the same as
        spec_fluxes.

    resampled_errs : numpy.ndarray
        Array of uncertainties associated with fluxes in
        res_fluxes. Only returned if spec_errs was specified.
    """

    old_spec_size = old_spec_wavs.shape[0]
    new_spec_size = new_spec_wavs.shape[0]

    # Arrays of left-hand sides and widths for the old bins
    spec_lhs = np.zeros(old_spec_size)
    spec_lhs[0] = old_spec_wavs[0] * 1.5 - old_spec_wavs[1] * 0.5
    spec_lhs[1:] = (old_spec_wavs[1:] + old_spec_wavs[:-1]) / 2

    spec_widths = np.zeros(old_spec_size)
    spec_widths[-1] = (old_spec_wavs[-1] - old_spec_wavs[-2])
    spec_widths[:-1] = spec_lhs[1:] - spec_lhs[:-1]

    # Arrays of left-hand sides and widths for the new bins
    filter_lhs = np.zeros(new_spec_size)
    filter_lhs[0] = new_spec_wavs[0] * 1.5 - new_spec_wavs[1] * 0.5
    filter_lhs[1:] = (new_spec_wavs[1:] + new_spec_wavs[:-1]) / 2

    filter_widths = np.zeros(new_spec_size)
    filter_widths[-1] = (new_spec_wavs[-1] - new_spec_wavs[-2])
    filter_widths[:-1] = filter_lhs[1:] - filter_lhs[:-1]

    # Arrays of the right-hand sides
    spec_rhs = spec_lhs + spec_widths[0]
    filter_rhs = filter_lhs + filter_widths[0]

    # New wavelength range is completely outside of the old range
    if filter_lhs[0] > spec_rhs[-1] or filter_rhs[-1] < spec_lhs[0]:
        raise ValueError(
            "spectres: New wavelength range is completely outside "
            "of the range of the old wavelength range.")

    # Part of the new wavelength range is outside of the old range
    # Only check and throw warning if silence is not set to True
    if not silence:
        if filter_lhs[0] < spec_lhs[0] or filter_rhs[-1] > spec_rhs[-1]:
            warnings.warn("spectres: Part of the new wavelength range "
                          "is outside the range of the input data, they are "
                          "filled with zeros.")

    # Generate output arrays to be populated
    res_fluxes = np.zeros(spec_fluxes[..., 0].shape + new_spec_wavs.shape)

    if spec_errs is not None:
        if spec_errs.shape != spec_fluxes.shape:
            raise ValueError("If specified, spec_errs must be the same shape "
                             "as spec_fluxes.")
        else:
            res_fluxerrs = np.copy(res_fluxes)

    # Find the first bin of the new spectrum
    if filter_lhs[0] < spec_lhs[0]:
        start_idx = np.amax(np.where(filter_lhs < spec_lhs[0]))
    else:
        start_idx = 0

    # Find the last bin of the new spectrum
    if filter_lhs[-1] > spec_rhs[-1]:
        stop_idx = np.amin(np.where(filter_lhs > spec_rhs[-1]))
    else:
        stop_idx = new_spec_size

    # Pad the new flux with the fill_value provided
    res_fluxes[:start_idx] += fill_value
    res_fluxes[stop_idx:] += fill_value

    # Indices of the old bins
    start = 0
    stop = 0

    # Calculate new flux and uncertainty values, loop over new bins
    for j in range(start_idx, stop_idx):

        # Find first old bin which is partially covered by the new bin
        while spec_rhs[start] < filter_lhs[j]:
            start += 1

        # Find last old bin which is partially covered by the new bin
        while ((spec_rhs[stop] < filter_rhs[j]) & (stop < old_spec_size - 1)):
            stop += 1

        # If new bin is fully within one old bin these are the same
        if stop == start:

            res_fluxes[..., j] = spec_fluxes[..., start]
            if spec_errs is not None:
                res_fluxerrs[..., j] = spec_errs[..., start]

        # Otherwise multiply the first and last old bin widths by P_ij
        else:

            start_factor = ((spec_rhs[start] - filter_lhs[j]) /
                            (spec_rhs[start] - spec_lhs[start]))

            end_factor = ((filter_rhs[j] - spec_lhs[stop]) /
                          (spec_rhs[stop] - spec_lhs[stop]))

            # Populate res_fluxes spectrum and uncertainty arrays
            spec_widths_ranged = spec_widths[start:stop + 1]
            spec_widths_ranged[0] *= start_factor
            spec_widths_ranged[-1] *= end_factor
            f_widths = spec_widths_ranged * spec_fluxes[..., start:stop + 1]
            res_fluxes[..., j] = np.sum(f_widths, axis=-1)
            res_fluxes[..., j] /= np.sum(spec_widths_ranged)

            if spec_errs is not None:
                e_wid = spec_widths_ranged * spec_errs[..., start:stop + 1]

                res_fluxerrs[..., j] = np.sqrt(np.sum(e_wid**2, axis=-1))
                res_fluxerrs[..., j] /= np.sum(spec_widths_ranged)

    # If errors were supplied return the res_fluxes spectrum and error arrays
    if spec_errs is not None:
        return res_fluxes, res_fluxerrs

    # Otherwise just return the res_fluxes spectrum array
    else:
        return res_fluxes
