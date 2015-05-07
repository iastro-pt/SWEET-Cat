#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
import numpy as np
import argparse
from astropy import coordinates as coord
from astropy import units as u
import csv
import pandas as pd
import operator
from clint.textui import puts, colored

from ParallaxSpec import parallax


def torres(name, teff=False, logg=False, feh=False):
    """
    Calculates the mass and error from Torres. See source for more information
    """
    from TorresMass import massTorres

    T, Terr = teff
    L, Lerr = logg
    F, Ferr = feh
    try:
        Terr, Lerr, Ferr = float(Terr), float(Lerr), float(Ferr)
        T, L, F = float(T), float(L), float(F)
    except ValueError:
        puts(colored.red('No mass derived for this star...'))
        return 'NULL', 'NULL'

    M, Merr = massTorres(T, Terr, L, Lerr, F, Ferr)
    puts(colored.green('Done'))
    return round(M, 2), round(Merr, 2)


def variable_assignment(digits=2):
    try:
        x = round(input('> '), digits)
    except SyntaxError, e:
        x = 'NULL'
    return x


def _parser():
    parser = argparse.ArgumentParser(description='Create a line with a new'
                                                 'update for SWEET-Cat'
                                                 'formatted for the webpage')
    parser.add_argument('input',
                        help='Name of the star as found in the exoplanet.eu'
                             'catalog')
    parser.add_argument('-o', '--output', default='newHost.rdb',
                        help='The output file. Note that only one line will be'
                             'written to this file, and it will be'
                             'overwritten')
    return parser.parse_args()


def _strip(text):
    """
    The star_name sometimes contains extra whitespaces. This function is parsed
    to pd.read_csv, and strip it from that column.
    """
    try:
        return text.strip()
    except AttributeError:
        return text


if __name__ == '__main__':
    args = _parser()
    input_, output = args.input.strip(), args.output
    # Read the data from exoplanet.eu
    fields = ['star_name', 'ra', 'dec', 'mag_v', 'star_metallicity',
              'star_teff']
    SC = pd.read_csv('exo.csv', skipinitialspace=True, usecols=fields,
                     converters={'star_name': _strip})

    SC = SC[SC['star_name'] == input_]

    try:
        name = SC.star_name.values[0]
    except IndexError, e:
        puts(colored.red(input_) + ' where not found. Try another star or add'
             'manually.')
        raise SystemExit

    # Get RA and dec
    ra, dec = float(SC.ra.values[0]), float(SC.dec.values[0])
    c = coord.SkyCoord(ra, dec, unit=(u.degree, u.degree), frame='icrs')
    RA = list(c.ra.hms)
    RA[0] = str(int(RA[0])).zfill(2)
    RA[1] = str(int(RA[1])).zfill(2)
    RA[2] = str(round(RA[2], 2)).zfill(4)
    if len(RA[2]) == 4:
        RA[2] += '0'
    RA = "{0} {1} {2}".format(*RA)

    DEC = list(c.dec.dms)
    DEC[0] = str(int(DEC[0])).zfill(2)
    DEC[1] = str(int(DEC[1])).zfill(2)
    DEC[2] = str(round(DEC[2], 2)).zfill(4)
    if len(DEC[2]) == 4:
        DEC[2] += '0'
    DEC = "+{0} {1} {2}".format(*DEC)

    # Here comes the user interface part...
    puts(colored.blue('\nStandard parameters\n\n'))

    # The HD number
    puts('The '+colored.yellow('HD number'))
    HD = raw_input('> ')
    if HD == '':
        HD = 'NULL'

    # The V magnitude
    V_exo = SC.mag_v.values[0]
    if np.isnan(V_exo):
        puts('The ' + colored.yellow('V magnitude'))
        V = variable_assignment(2)
        puts('The error on ' + colored.yellow('V magnitude'))
        Verr = variable_assignment(2)
    else:
        V = round(float(V_exo), 2)
        puts('The error on ' + colored.yellow('V magnitude'))
        Verr = variable_assignment(2)

    # The metallicity
    FeH_exo = SC.star_metallicity.values[0]
    if np.isnan(FeH_exo):
        puts('The ' + colored.yellow('[Fe/H]'))
        FeH = variable_assignment(2)
        puts('The error on ' + colored.yellow('[Fe/H]'))
        Ferr = variable_assignment(2)
    else:
        FeH = round(float(FeH_exo), 2)
        puts('The error on ' + colored.yellow('[Fe/H]'))
        Ferr = variable_assignment(2)

    # The effective temperature
    Teff_exo = SC.star_teff.values[0]
    if np.isnan(Teff_exo):
        puts('The ' + colored.yellow('Teff'))
        Teff = variable_assignment(0)
        puts('The error on ' + colored.yellow('Teff'))
        Tefferr = int(variable_assignment(0))
    else:
        Teff = int(Teff_exo)
        puts('The error on ' + colored.yellow('Teff'))
        Tefferr = int(variable_assignment(0))

    # The log g
    puts('The ' + colored.yellow('logg'))
    logg = variable_assignment(2)
    puts('The error on ' + colored.yellow('logg'))
    loggerr = variable_assignment(2)

    # The mass
    puts(colored.magenta('\nCalculating the mass...'))
    M, Merr = torres(name, [Teff, Tefferr], [logg, loggerr], feh=[FeH, Ferr])

    # The parallax
    puts('\nIs the ' + colored.yellow('parallax')+' given from SIMBAD?')
    par = raw_input('(y/n) > ')
    if par.lower() == 'y' or par.lower() == 'yes':
        puts('The '+colored.yellow('parallax'))
        p = input('> ')
        puts('The error on '+colored.yellow('parallax'))
        perr = input('> ')
        pflag = 'Simbad'
    elif par.lower() == 'n' or par.lower() == 'no' or par == '':
        puts(colored.magenta('\nCalculating the parallax...'))
        try:
            p = round(parallax(Teff, logg, V, M), 2)
            puts(colored.green('Done'))
        except TypeError:
            puts(colored.red('Could not calculate the parallax...'))
            p = 'NULL'
        perr = 'NULL'
        pflag = 'Spec'
    else:
        p = 'NULL'
        perr = 'NULL'
        pflag = 'NULL'
        puts(colored.red('Parallax, the error, and the flag all set to NULL'))

    # The microturbulence number
    puts('The '+colored.yellow('microturbulence'))
    vt = variable_assignment(2)
    puts('The error on '+colored.yellow('microturbulence'))
    vterr = variable_assignment(2)

    # Author and link to ADS
    puts('Who is the '+colored.yellow('author?'))
    author = raw_input('> ')
    puts('Link to article ('+colored.yellow('ADS')+')')
    link = raw_input('> ')

    # Source flag
    puts(colored.yellow('Source flag'))
    source = raw_input('(0/1) > ')
    if source == '':
        source = 'NULL'

    # Last update
    puts(colored.yellow('Last update'))
    update = raw_input('(1989-09-13) > ')

    # Comments
    puts('Any '+colored.yellow('comments'))
    puts('E.g. if we have a M dwarf...')
    comment = raw_input('> ')
    if comment == '':
        comment = 'NULL'

    params = [name, HD, RA, DEC, V, Verr, p, perr, pflag, Teff, Tefferr,
              logg, loggerr, 'NULL', 'NULL', vt, vterr, FeH, Ferr, M, Merr,
              author, link, source, update, comment]
    params = map(str, params)

    with open(output, 'wb') as f:
        f.write('\t'.join(params) + '\tNULL\n')
