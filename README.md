# SpectRes

SpectRes is a Python function which efficiently resamples spectra and their associated uncertainties onto an arbitrary wavelength grid. The function works with any grid of wavelength values, including non-uniform sampling, and preserves the integrated flux. 

This may be of use for binning data to increase the signal to noise ratio, obtaining synthetic photometry, or resampling model spectra to match the sampling of observed data for spectral energy distribution fitting. 

The function takes three arguments, firstly __**spec_wavs**__, an array of wavelengths corresponding to the current sampling of the spectrum. Secondly __**spec_fluxes**__, a 1D or 2D array of flux values with the first axis running over wavelength and the second (if present) over different spectra to be resampled. The third argument, __**resampling**__, is the desired sampling. A keyword argument, __**spec_errs**__ may also be passed of the same shape as __**spec_fluxes**__ containing the uncertainty on each flux value. 

The function returns an array, __**resampled**__ containing the resampled spectrum or spectra, with first dimension the same length as __**resampling**__ and second dimension (if present) the same length as the second dimension of __**spec_fluxes**__. If the keyword argument __**spec_errs**__ is passed, a second array, __**resampled_errs**__ of the same shape, containing the resampled error spectrum or spectra is also returned.
