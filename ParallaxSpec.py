#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def bolcor(teff):
    """
    Calculate the bolometric correction, given the temperature
    """
    lteff = np.log10(teff)
    bcflow=np.nan

    if lteff < 3.7:
        bcflow = -0.190537291496456e+05 + 0.155144866764412e+05 * lteff -\
            0.421278819301717e+04 * (lteff * lteff) +\
            0.381476328422343e+03 * (lteff*lteff*lteff)

    if lteff >= 3.7 and lteff < 3.9:
        bcflow = -0.370510203809015e+05 + 0.385672629965804e+05 * lteff -\
            0.150651486316025e+05 * (lteff * lteff) +\
            0.261724637119416e+04 * (lteff*lteff*lteff) -\
            0.170623810323864e+03 * (lteff * lteff * lteff * lteff)
    else:
        bcflow = -0.118115450538963e+06 + 0.137145973583929e+06 * lteff -\
            0.636233812100225e+05 * (lteff * lteff) +\
            0.147412923562646e+05 * (lteff * lteff * lteff) -\
            0.170587278406872e+04 * (lteff * lteff * lteff * lteff) +\
            0.788731721804990e+02 * (lteff * lteff * lteff * lteff * lteff)
    
    return bcflow


def parallax(teff,eteff, logg,elogg,vmag,evmag,   mass,emass,  Av,eAv):
    """
    Calculate the parallax, given the mass Santos 2004
    """
    if eteff=='NULL' or elogg=='NULL' or evmag=='NULL' or  emass=='NULL' or  eAv=='NULL':
        bcflow = bolcor(teff)
        return 10.**((logg - 4.44 - np.log(mass)/np.log(10.) - 4.*np.log(teff) / np.log(10.) + \
            4.*np.log(5777.)/np.log(10.) - 0.4*(vmag + bcflow-Av) - 0.11) * 0.5) * 1000, np.nan

    ntrials = 10000
    randomteff = float(teff) + float(eteff)*np.random.randn(ntrials)
    randomlogg = float(logg) + float(elogg)*np.random.randn(ntrials)
    randommass = abs(mass + emass*np.random.randn(ntrials))
    randomvmag = float(vmag) + float(evmag)*np.random.randn(ntrials)
    randomAv = Av + eAv*np.random.randn(ntrials)
    par=[]
    for i in range(ntrials):
        bcflow = bolcor(randomteff[i])
        par.append(10.**((randomlogg[i] - 4.44 - np.log(randommass[i])/np.log(10.) \
        - 4.*np.log(randomteff[i]) / np.log(10.) + 4.*np.log(5777.)/np.log(10.) - 0.4*\
        (randomvmag[i] + bcflow-randomAv[i]) - 0.11) * 0.5) * 1000)
    
    sig=np.sqrt(np.sum((par - np.mean(par))**2) / (ntrials - 1))
    return np.mean(par),sig    
