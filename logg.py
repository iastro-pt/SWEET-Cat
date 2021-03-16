#!/usr/bin/python

import argparse
from astropy import constants as c
from numpy import log10
import numpy as np

def _parse():
    '''Calculate the surface gravity from M and R'''
    p = argparse.ArgumentParser(description='Calculate logg from solar M and R')
    p.add_argument('M', help='Mass in solar units', type=float)
    p.add_argument('R', help='Radius in solar units', type=float)
    return p.parse_args()


def logg(M, R):
    """Mass and radius in solar units"""
    G = c.G.value * 1e3
    M *= c.M_sun.value * 1e3
    R *= c.R_sun.value * 1e2
    return log10(G*M/(R**2))

def BC_g_Andrae(Teff, Tsun = 5777):
    a  = [6e-2    , 6.731e-5, -6.647e-8,  2.859e-11, -7.197e-15]
    sa = [2.634e-2, 2.438e-5, -1.129e-9, -6.722e-12,  1.635e-15]
    if Teff < 4000:
        a  = [ 1.749,  1.977e-3, 3.737e-7, -8.966e-11, -4.183e-14]
        sa = [-2.487, -1.8762-3, 2.128e-7,  3.807e-10,  6-570e-14]

    bcg  = np.sum([a[i]*(Teff-Tsun)**i for i in range(len(a))])
    sbcg = np.sum([sa[i]*(Teff-Tsun)**i for i in range(len(sa))])
    return bcg, sbcg

def logg_gaia(mass, teff, gaia_gmag, gaia_paralax, Ag = 0, teff_sun=5777):

    MG = gaia_gmag + 5 - 5*log10(1000./gaia_paralax) - Ag
    bcg,sbcg = BC_g_Andrae(teff)
    mbol_star = MG + bcg
    mbol_sun = 4.74  #IAU Resolution 2015
    logg_sun =  logg(1,1)
    logg_gaia = log10(mass) + 4*log10(teff)-4*log10(teff_sun) + 0.4*mbol_star - 0.4*mbol_sun + logg_sun
    return logg_gaia

def main():
    args = _parse()
    print(('logg: %.2f' % logg(args.M, args.R)))


if __name__ == '__main__':
    main()
