import numpy as np
from timeit import timeit, Timer

from spectral_resampling_numba import spectres as spnumba
from spectral_resampling import spectres as sp

def call_spectres():
    inwaves = np.linspace(0.5, 15, size)
    outwaves = np.linspace(1, 10, 200)
    influxes = np.ones_like(inwaves)
    r = sp(outwaves, inwaves, influxes)
    return r

def call_spectres_err():
    inwaves = np.linspace(0.5, 15, size)
    outwaves = np.linspace(1, 10, 200)
    influxes = np.ones_like(inwaves)
    inerrs = 0.1*influxes
    r = sp(outwaves, inwaves, influxes, spec_errs = inerrs)
    return r


def call_spectresn():
    inwaves = np.linspace(0.5, 15, size)
    outwaves = np.linspace(1, 10, 200)
    influxes = np.ones_like(inwaves)
    r = spnumba(outwaves, inwaves, influxes)
    return r

def call_spectresn_err():
    inwaves = np.linspace(0.5, 15, size)
    outwaves = np.linspace(1, 10, 200)
    influxes = np.ones_like(inwaves)
    inerrs = 0.1*influxes
    r = spnumba(outwaves, inwaves, influxes, spec_errs = inerrs)
    return r


if __name__=="__main__":
    insizes = [1000,10000, 100000]
    repeats = 1000
    size = insizes[0]
    a = call_spectres_err()

    b = call_spectresn_err()

    for size in insizes:
        
        print("Running spectres with ",size," input wavelengths for ",repeats," times")
        timer = Timer(stmt='call_spectres()', setup = 'from spectral_resampling import spectres as sp; from __main__ import call_spectres') 
        time = timer.timeit(repeats)
        print("Total runtime = ",time,"s, time per call = ",time/repeats,"s")

        call_spectresn() #Call the compiled version once with the same argument types so that it's ready for the speed test - this could really be done outside the loop
        print("Running compiled version of spectres with ",size," input wavelengths for ",repeats," times")
        timer = Timer(stmt='call_spectresn()', setup = 'from spectral_resampling import spectres as sp; from __main__ import call_spectresn')
        time_c = timer.timeit(repeats)

        print("Total runtime = ",time_c,"s, time per call = ",time_c/repeats,"s")
        print("Speedup thanks to Numba = ",time/time_c)

        
        

    
