SpectRes: Simple Spectral Resampling
====================================

SpectRes is a Python function which efficiently resamples spectra and their associated uncertainties onto an arbitrary wavelength grid. The function works with any grid of wavelength values, including non-uniform sampling, and preserves the integrated flux. This may be of use for binning data to increase the signal to noise ratio, obtaining synthetic photometry, or resampling model spectra to match the sampling of observational data. 


Installation
------------

SpectRes can be installed using pip 

.. code::

	pip install spectres

	
Documentation
-------------

The code's documentation is available at `spectres.readthedocs.io <https://spectres.readthedocs.io>`_.  Take a look in the examples folder of this repository for some more guidance.

For an explanation of how the code works take a look at `ArXiv:1705.05165 <https://arxiv.org/abs/1705.05165>`_. Please consider citing this publication if you use SpectRes in your research.