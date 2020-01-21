#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Example:
# python SC_exoplanet_v2.py teq0 sma -z metal -lx -ly -o star teq0 sma metal -t aaa.rdb
#

# My imports

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
try:
    import seaborn as sns
    sns.set_style('dark')
    sns.set_context('talk', font_scale=1.2)
    color = sns.color_palette()
except ImportError:
    print('Install seaborn for better plots (optional): pip install seaborn')
    color = 'b,g,r,m,y,k'.split(',')
try:
    import pandas as pd
except ImportError:
    raise ImportError('Install pandas: pip install pandas')
from PyAstronomy import pyasl
import argparse


def radTorres(teff, erteff, logg, erlogg, feh, erfeh):
    ntrials = 100
    randomteff = teff + erteff*np.random.randn(ntrials)
    randomlogg = logg + erlogg*np.random.randn(ntrials)
    randomfeh = feh + erfeh*np.random.randn(ntrials)

    # Parameters for the Torres calibration:
    b1, b2, b3 = 2.4427, 0.6679, 0.1771
    b4, b5, b6 = 0.705, -0.21415, 0.02306
    b7 = 0.04173

    logR = np.zeros(ntrials)
    for i in range(ntrials):
        X = np.log10(randomteff[i]) - 4.1
        logR[i] = b1 + b2*X + b3*X**2 + b4*X**3 + b5*randomlogg[i]**2 + b6*randomlogg[i]**3 + b7*randomfeh[i]

    meanRadlog = np.mean(logR)
    sigRadlog = np.sqrt(np.sum((logR-meanRadlog)**2))/(ntrials-1)
    sigRadlogTot = np.sqrt(0.014**2 + sigRadlog**2)

    meanRad = 10**meanRadlog
    sigRad = 10**(meanRadlog + sigRadlogTot) - meanRad
    return meanRad, sigRad


def _parser():
    parser = argparse.ArgumentParser(description='Preprocess the results')
    sc = pyasl.SWEETCat()
    sc = sc.data
    eu = pyasl.ExoplanetEU()
    eu = eu.getAllData()
    eu = pd.DataFrame(eu)
    p = list(sc.columns) + list(eu.columns) + ['teq0', 'mass', 'masserr', 'radius', 'radiuserr', 'lum']
    parser.add_argument('x', choices=p)
    parser.add_argument('y', choices=p)
    parser.add_argument('-z', help='Color scale', choices=p, default=None)
    parser.add_argument('-i', '--input', help='File name of result file', default='results.csv')
    parser.add_argument('-c', '--convergence', help='Only plot converged results', default=True, action='store_false')
    parser.add_argument('-ix', help='Inverse x axis', default=False, action='store_true')
    parser.add_argument('-iy', help='Inverse y axis', default=False, action='store_true')
    parser.add_argument('-lx', help='Logarithmic x axis', default=False, action='store_true')
    parser.add_argument('-ly', help='Logarithmic y axis', default=False, action='store_true')
    parser.add_argument('-o', '--output', help='Save a list of arguments to exoplanets.csv', nargs='+')
    parser.add_argument('-t', '--table', help='Table to intercept. Column name must be "star"')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = _parser()
    # Prepare the SWEET-Cat data
    print('Downloading the data from SWEET-Cat...')
    sc = pyasl.SWEETCat()
    sc.downloadData()
    sc = sc.data

    newname = [x.lower().replace(' ', '') for x in sc.star]  # Put name lower and remove all spaces
    newname = list(map(str.strip, newname))  # Remove all newline or tab characters
    sc['nameNew'] = newname

    # Prepare the exoplanetEU data
    print('Downloading the data from exoplanetEU...')
    eu = pyasl.ExoplanetEU()
    eu = eu.getAllData()
    eu = pd.DataFrame(eu)  # Convert the structure to a DataFrame

    newstname = [x.lower().replace(' ', '') for x in eu.stName]  # Put name lower and remove all spaces
    newstname = list(map(str.strip, newstname))  # Remove all newline or tab characters
    eu['stNameNew'] = newstname

    # Merge the two based on the stellar name
    print(1)
    df = pd.merge(left=sc, right=eu, left_on='nameNew', right_on='stNameNew')
    print(2)
    df.rename(columns={'ra_x': 'ra', 'dec_x': 'dec'}, inplace=True)

    # Calculate radius and luminosity
    rr = ['radius', 'radiuserr', 'teq0']
    if (args.x in rr) or (args.y in rr) or (args.z in rr):
        params = list(zip(df.teff, df.erteff, df.logg, df.erlogg, df.metal, df.ermetal))
        print(3)
        r = [radTorres(t, et, l, el, f, ef) for t, et, l, el, f, ef in params]
        # Insert the radius, and luminosity in the table
        print(4)
        df['radius'] = pd.Series(np.asarray(r)[:, 0])
        df['radiuserr'] = pd.Series(np.asarray(r)[:, 1])
        # Compute Teq0:
        df['teq0'] = df.teff*((df.radius*700000.)/(2.*df.sma*150000000.))**(0.5)

    df['lum'] = (df.teff/5777.)**4 * df.mass

    # Intersect the table df with the list of stars in file aaa.rdb
    if args.table:
        tt = pd.read_csv(args.table)
        df2 = pd.merge(left=tt, right=df, left_on='planet', right_on='plName', how='inner')
        df = df2

    plt.figure()
    if args.z:
        z = df[args.z].values
        color = df[args.z].values
        z[np.isnan(z)] = min(z[~np.isnan(z)])
        size = (z-z.min())/(z.max()-z.min())*100
        size[np.argmin(size)] = 10  # Be sure to see the "smallest" point
        plt.scatter(df[args.x], df[args.y], c=color, s=size, cmap=cm.seismic)
    else:
        plt.scatter(df[args.x], df[args.y], c=color[0], s=40)

    labels = {'teff': r'$T_\mathrm{eff}$ [K]',
              'tefferr': r'$\sigma T_\mathrm{eff}$ [K]',
              'stTeff': r'$T_\mathrm{eff}$ [K]',
              'logg': r'$\log(g)$ [cgs]',
              'logglc': r'$\log(g)$ [cgs]',
              'loggerr': r'\sigma $\log(g)$ [cgs]',
              'erlogglc': r'\sigma $\log(g)$ [cgs]',
              'metal': '[Fe/H]',
              'ermetal': r'$\sigma$ [Fe/H]',
              'vt': r'$\xi_\mathrm{micro}$ [km/s]',
              'vterr': r'$\sigma\xi_\mathrm{micro}$ [km/s]',
              'lum': r'$L_\odot$',
              'mass': r'$M_\odot$',
              'ermass': r'$\sigma M_\odot$',
              'stMass': r'$M_\odot$',
              'stRadius': r'$R_\odot$',
              'sma': 'Semi major axis [AU]',
              'vmag': 'V magnitude',
              'ervmag': 'Error on V magnitude',
              'par': 'Parallax ["]',
              'erpar': 'Error on parallax ["]',
              'plMass': r'Planet mass [$M_\mathrm{Jup}$]',
              'plRadius': r'Planet radius [$R_\mathrm{Jup}$]',
              'period': 'Period [day]',
              'eccentricity': 'Eccentricity',
              'inclination': 'Inclination [degree]',
              'angDistance': 'Angular distance [arcsec]',
              'discovered': 'Discovered [year]',
              'omega': r'$\Omega$ [degree]',
              'tperi': 'Time of periastron [JD]',
              'mag_v': 'V mag',
              'mag_i': 'I mag',
              'mag_j': 'J mag',
              'mag_h': 'H mag',
              'mag_k': 'K mag',
              'dist': 'Distance [pc]',
              'mh': 'Metallicity',
              'stAge': 'Stellar age [Gyr]'}

    try:
        plt.xlabel(labels[args.x])
    except KeyError:
        plt.xlabel(args.x)
    try:
        plt.ylabel(labels[args.y])
    except KeyError:
        plt.ylabel(args.y)

    if args.z:
        cbar = plt.colorbar()
        try:
            cbar.set_label(labels[args.z])
        except KeyError:
            cbar.set_label(args.z)

    if args.lx:
        plt.xscale('log')
    if args.ly:
        plt.yscale('log')
    if args.ix:
        plt.xlim(plt.xlim()[::-1])
    if args.iy:
        plt.ylim(plt.ylim()[::-1])

    plt.tight_layout()
    plt.grid(True)
    plt.show()

    if args.output:
        try:
            dfout = df[args.output]
        except KeyError:
            for key in args.output:
                try:
                    _ = df[key]
                except KeyError:
                    print('Invalid name: %s' % key)
            raise SystemExit()

        dfout.to_csv('exoplanets.csv', sep='\t', index=False, na_rep='...')
        print('Saved result in exoplanet.csv')
