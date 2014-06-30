#!/usr/bin/env python
# encoding: utf-8

import sys
import numpy as np


inputPAR = sys.argv[1]


def bolcor(teff):
    """
    Calculate the bolometric correction, given the temperature
    """
    lteff = np.log10(teff)

    if lteff < 3.7:
        bcflow = -0.190537291496456e+05 + 0.155144866764412e+05 * lteff -\
            0.421278819301717e+04 * (lteff * lteff) +\
            0.381476328422343e+03 * (lteff*lteff*lteff)
        return bcflow

    if lteff >= 3.7 and lteff < 3.9:
        bcflow = -0.370510203809015e+05 + 0.385672629965804e+05 * lteff -\
            0.150651486316025e+05 * (lteff * lteff) +\
            0.261724637119416e+04 * (lteff*lteff*lteff) -\
            0.170623810323864e+03 * (lteff * lteff * lteff * lteff)
        return bcflow
    else:
        bcflow = -0.118115450538963e+06 + 0.137145973583929e+06 * lteff -\
            0.636233812100225e+05 * (lteff * lteff) +\
            0.147412923562646e+05 * (lteff * lteff * lteff) -\
            0.170587278406872e+04 * (lteff * lteff * lteff * lteff) +\
            0.788731721804990e+02 * (lteff * lteff * lteff * lteff * lteff)
        return bcflow


def parallax(teff, logg, vmag, mass, bcflow):
    """
    Calculate the parallax, given the mass
    """
    par = 10.**((logg - 4.44 - np.log(mass)/np.log(10.) - 4.*np.log(teff) /
                np.log(10.) + 4.*np.log(5777.)/np.log(10.) - 0.4*(vmag +
                                                                  bcflow) -
                0.11) * 0.5) * 1000

    return par


with open(inputPAR) as f:
    lines = f.readlines()

for line in lines:
    words = line.split('    ')

    name = words[0]
    teff = float(words[1])
    logg = float(words[2])
    mass = float(words[3])
    vmag = float(words[4])

    bcflow = bolcor(teff)
    par = parallax(teff, logg, vmag, mass, bcflow)
    print "Name        parallax"
    print "===================="
    print name+':  ', round(par, 2)
