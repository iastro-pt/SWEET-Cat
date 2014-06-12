#!/usr/bin/python
#
# Program to find the Stellar Mass through Torres
#
# The program needs one inputfile (rdb format):
# File format of stellar atmosphere parameters: Name Teff errTeff Logg errLogg
# Metal errMetal (...)
#
# Usage: python TorresMass.py file.rdb
#
# Program will print the output on the screen

import sys
import numpy as np

inputMOOG = sys.argv[1]


# Calculate a mass using the Torres calibration
def massTorres(teff, erteff, logg, erlogg, feh, erfeh):

    ntrials = 10000
    randomteff = teff + erteff * np.random.randn(ntrials)
    randomlogg = logg + erlogg * np.random.randn(ntrials)
    randomfeh = feh + erfeh * np.random.randn(ntrials)

    # Parameters for the Torres calibration:
    a1 = 1.5689
    a2 = 1.3787
    a3 = 0.4243
    a4 = 1.139
    a5 = -0.1425
    a6 = 0.01969
    a7 = 0.1010

    M = np.zeros(ntrials)
    logM = np.zeros(ntrials)

    for i in xrange(len(randomteff)):

        X = np.log10(randomteff[i]) - 4.1
        logMass = a1 + a2 * X + a3 * X * X + a4 * X * X * X + a5 *\
            randomlogg[i] * randomlogg[i] + a6 * randomlogg[i] *\
            randomlogg[i] * randomlogg[i] + a7 * randomfeh[i]
        logM[i] = logMass
        M[i] = 10 ** logMass

    meanMasslog = np.mean(logM)
    sigMasslog = np.sqrt(np.sum([(logMi - meanMasslog)**2 for logMi in logM]) /
                         (ntrials - 1))
    sigMasslogTot = np.sqrt(0.027*0.027 + sigMasslog*sigMasslog)

    meanMass = 10**meanMasslog
    sigMass = 10**(meanMasslog + sigMasslogTot) - meanMass

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

    for i in range(len(randomteff)):

        X = np.log10(randomteff[i]) - 4.1
        logRad = b1 + b2 * X + b3 * X * X + b4 * X * X * X + b5 *\
            randomlogg[i] * randomlogg[i] + b6 * randomlogg[i] *\
            randomlogg[i] * randomlogg[i] + b7 * randomfeh[i]
        logR[i] = logRad
        R[i] = 10 ** logRad

    meanRadlog = np.mean(logR)
    sigRadlog = np.sqrt(np.sum([(logRi-meanRadlog)**2 for logRi in logR]) /
                        (ntrials-1))
    sigRadlogTot = np.sqrt(0.014*0.014 + sigRadlog*sigRadlog)

    meanRad = 10**meanRadlog
    sigRad = 10**(meanRadlog + sigRadlogTot) - meanRad

    return meanRad, sigRad


######################################################################
# Main program
######################################################################


with open(inputMOOG) as f:

    f.readline()
    f.readline()
    lines = f.readlines()

p = [0.79102951, 0.0 - 0.57480278,  0.70062477]

print "Name,  mass,   errmass"
print "Name, radius, errRadius"
print "======================="
for line in lines:
    words = line.split('\t')
    name = words[0]

    teff = float(words[1])
    erteff = float(words[2])
    logg = float(words[3])
    erlogg = float(words[4])
    feh = float(words[5])
    erfeh = float(words[6])

    M = massTorres(teff, erteff, logg, erlogg, feh, erfeh)

    mass = str(round(p[0]*M[0]*M[0] + p[1]*M[0] + p[2], 2))
    ermass = str(round(M[1], 2))

    print name, mass, ermass

    rad = 1
    if rad:
        R = radTorres(teff, erteff, logg, erlogg, feh, erfeh)
        rad = str(round(R[0], 2))
        errad = str(round(R[1], 2))

        print name, rad, errad
