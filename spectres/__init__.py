try:
    from .spectral_resampling_numba import spectres
except ImportError:
    from .spectral_resampling import spectres
