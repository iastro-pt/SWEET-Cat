#!/usr/bin/python

import argparse
from astropy import constants as c
from numpy import log10


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


def main():
    args = _parse()
    print 'logg: %.2f' % logg(args.M, args.R)


if __name__ == '__main__':
    main()
