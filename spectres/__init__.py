from .spectral_resampling import spectres

try:
    import numba
    from .spectral_resampling_numba import spectres_numba

except ImportError:
    pass
