.. SpectRes documentation master file, created by
   sphinx-quickstart on Sat Mar 31 10:52:51 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

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

The code is developed at `github.com/ACCarnall/spectres <https://github.com/ACCarnall/spectres>`_. Take a look in the examples folder for some more guidance.

For an explanation of how the code works take a look at `ArXiv1705.05165 <https://arxiv.org/abs/1705.05165>`_. Please consider citing this publication if you use SpectRes in your research.

Numba compiled version
----------------------

With thanks to Peter Scicluna, SpectRes now comes with an optional Numba compiled version, which should speed the code up under a range of circumstances. You can try this out by calling the ``spectres.spectres_numba`` function in place of ``spectres.spectres``.

API Documentation
-----------------

.. autofunction:: spectres.spectres
