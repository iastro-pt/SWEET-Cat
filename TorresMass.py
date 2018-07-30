#!/usr/bin/python

import numpy as np
from astroquery.vizier import Vizier

def massTorres(teff, erteff, logg, erlogg, feh, erfeh):
    """ 
    Calculate stellar mass using the Torres et al. (2010) callibration.
    Parameters
    ----------
    teff, erteff : floats
        Effective temperature and associated uncertainty.
    logg, erlogg : floats
        Surface gravity and associated uncertainty.
    feh, erfeh : floats
        Metallicity [Fe/H] and associated uncertainty.
   
    Returns
    -------
    meanMass, sigMass : floats
        Estimate for the stellar mass and associated uncertainty.
    """

    # Number of Monte Carlo trials for the uncertainty calculation.
    ntrials=10000

    randomteff = teff + erteff * np.random.randn(ntrials)
    randomlogg = logg + erlogg * np.random.randn(ntrials)
    randomfeh = feh + erfeh * np.random.randn(ntrials)

    # Parameters for the Torres calibration
    a1 = 1.5689
    a2 = 1.3787
    a3 = 0.4243
    a4 = 1.139
    a5 = -0.1425
    a6 = 0.01969
    a7 = 0.1010

    X = np.log10(randomteff) - 4.1
    
    logMass = a1 + a2*X + a3*X**2 + a4*X**3 + a5*randomlogg**2 + a6*randomlogg**3 + a7*randomfeh
    meanlogMass = np.mean(logMass)
    siglogMass = np.sum((logMass - meanlogMass)**2) / (ntrials - 1)
    # Add (quadratically) the intrinsic error of the calibration (0.027 in log mass).
    siglogMass = np.sqrt(0.027**2 + siglogMass)
    meanMass = 10**meanlogMass
    sigMass = 10**(meanlogMass + siglogMass) - meanMass
    
    # Correct the mass for the offset relative to isochrone-derived masses.
    if .7<=meanMass<=1.3:
        # correction comes from Santos+(2013), the SWEET-Cat paper
        randomMass = meanMass + sigMass * np.random.randn(ntrials)
        corrected_Mass = 0.791 * randomMass**2 - 0.575 * randomMass + 0.701
        meanMassCor = np.mean(corrected_Mass)
        sigMassCor = np.sqrt(np.sum((corrected_Mass - meanMassCor)**2) / (ntrials - 1))
        return meanMassCor, sigMassCor

    return meanMass, sigMass


def radTorres(teff, erteff, logg, erlogg, feh, erfeh):
    ntrials = 10000
    randomteff = teff + erteff*np.random.randn(ntrials)
    randomlogg = logg + erlogg*np.random.randn(ntrials)
    randomfeh = feh + erfeh*np.random.randn(ntrials)

    # Parameters for the Torres calibration:
    b1 = 2.4427
    b2 = 0.6679
    b3 = 0.1771
    b4 = 0.705
    b5 = -0.21415
    b6 = 0.02306
    b7 = 0.04173

    R = np.zeros(ntrials)
    logR = np.zeros(ntrials)

    X = np.log10(randomteff) - 4.1
    logRad = b1 + b2 * X + b3 * X * X + b4 * X * X * X + b5 * randomlogg * randomlogg \
    + b6 * randomlogg * randomlogg * randomlogg + b7 * randomfeh
    R = 10 ** logRad

    meanRadlog = np.mean(logRad)
    sigRadlog = np.sqrt(np.sum((logRad-meanRadlog)**2) / (ntrials-1))
    sigRadlog = np.sqrt(0.014**2 + sigRadlog**2)

    meanRad = 10**meanRadlog
    sigRad = 10**(meanRadlog + sigRadlog) - meanRad

    return meanRad, sigRad
