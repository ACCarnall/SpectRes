# SpectRes

SpectRes is a Python function which efficiently resamples spectra and their associated uncertainties onto an arbitrary wavelength grid. The function works with any grid of wavelength values, including non-uniform sampling, and preserves the integrated flux. 

This may be of use for binning data to increase the signal to noise ratio, obtaining synthetic photometry, or resampling model spectra to match the sampling of observed data for spectral energy distribution fitting. 

An article explaining the method the code employs will be published at www.ArXiv.org in the coming days. If you make use of the code in your research, please cite this article in any publications.

## Installation

The code is essentially pretty simple and so is currently formatted as a single Python function. If you want the function to be available to your other Python scripts, either copy the _spectres.py_ file into the directory containing the scripts, or you can add the directory you clone the repository into to your **PYTHONPATH** variable.

## Using The Code

The function takes three arguments, firstly _**spec_wavs**_, an array of wavelengths corresponding to the current sampling of the spectrum. Secondly _**spec_fluxes**_, a 1D or 2D array of flux values with the first axis running over wavelength and the second (if present) over different spectra to be resampled. The third argument, _**resampling**_, is the desired sampling. A keyword argument, _**spec_errs**_ may also be passed of the same shape as _**spec_fluxes**_ containing the uncertainty on each flux value. 

The function returns an array, _**resampled**_ containing the resampled spectrum or spectra, with first dimension the same length as _**resampling**_ and second dimension (if present) the same length as the second dimension of _**spec_fluxes**_. If the keyword argument _**spec_errs**_ is passed, a second array, _**resampled_errs**_ of the same shape, containing the resampled error spectrum or spectra is also returned.

## Example Files

Two examples are provided, the first resamples the spectrum of a high redshift quasar from Carnall et al. (2015) onto a coarser wavelength grid (bins the data) in order to improve the signal to noise ratio. The second resamples a whole grid of BC03 models (available at http://www.bruzual.org/bc03) onto a new wavelength grid.
