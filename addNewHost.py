#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
import numpy as np
import argparse
from astropy import coordinates as coord
from astropy import units as u
import csv
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a line with a new\
    update for SWEET-Cat formatted for the webpage')

    parser.add_argument('-i',
                        required=True,
                        help='Name of the star as found in the exoplanet.eu\
                              catalog')
    parser.add_argument('-o',
                        default='newHost.rdb',
                        help='The output file. Note that only one line will be\
                              written to this file, and it will be\
                              overwritten')

    args = parser.parse_args()

    SC = csv.reader(open('exo.csv'), delimiter=',')
    SC = sorted(SC, key=operator.itemgetter(24), reverse=True)
    n = len(SC)
    starID = ['']*n
    for i in xrange(n):
        starID[i] = SC[i][48]
    starID = np.array(starID)
    idx = np.where(args.i == starID)[0]

    if len(idx) > 1:  # For multiple systems
        idx = np.array(idx[0])

    t = SC[idx]
    name = t[48]
    ra, dec = float(t[49]), float(t[50])
    c = coord.ICRS(ra, dec, unit=(u.degree, u.degree))
    RA = list(c.ra.hms)
    RA[0] = str(int(RA[0])).zfill(2)
    RA[1] = str(int(RA[1])).zfill(2)
    RA[2] = str(round(RA[2], 2)).zfill(4)
    if len(RA[2]) == 4:
        RA[2] += '0'
    RA = RA[0] + " " + RA[1] + " " + RA[2]

    DEC = list(c.dec.dms)
    DEC[0] = str(int(DEC[0])).zfill(2)
    DEC[1] = str(int(DEC[1])).zfill(2)
    DEC[2] = str(round(DEC[2], 2)).zfill(4)
    if len(DEC[2]) == 4:
        DEC[2] += '0'
    DEC = '+' + DEC[0] + " " + DEC[1] + " " + DEC[2]

    # Here comes the user interface part...
    puts(colored.blue('\nStandard parameters\n\n'))

    # The HD number
    puts('The '+colored.yellow('HD number'))
    HD = raw_input('> ')
    if HD == '':
        HD = 'NULL'

    # The V magnitude
    if t[51] == '':
        puts('The ' + colored.yellow('V magnitude'))
        V = round(input('> '), 2)
        puts('The error on' + colored.yellow('V magnitude'))
        Verr = raw_input('> ')
        if Verr == '':
            Verr = 'NULL'
    else:
        V = round(float(t[51]), 2)
        puts('The error on ' + colored.yellow('V magnitude'))
        Verr = raw_input('> ')
        if Verr == '':
            Verr = 'NULL'


# The metallicity
    if t[57] == '':
        puts('The ' + colored.yellow('[Fe/H]'))
        FeH = round(input('> '), 2)
        puts('The error on ' + colored.yellow('[Fe/H]'))
        Ferr = raw_input('> ')
        if Ferr == '':
            Ferr = 'NULL'
    else:
        FeH = round(float(t[57]), 2)
        puts('The error on ' + colored.yellow('[Fe/H]'))
        Ferr = raw_input('> ')
        if Ferr == '':
            Ferr = 'NULL'


# The effective temperature
    if t[62] == '':
        puts('The ' + colored.yellow('Teff'))
        Teff = input('> ')
        puts('The error on ' + colored.yellow('Teff'))
        Tefferr = raw_input('> ')
        if Tefferr == '':
            Tefferr = 'NULL'
    else:
        Teff = int(float(t[62]))
        puts('The error on ' + colored.yellow('Teff'))
        Tefferr = raw_input('> ')
        if Tefferr == '':
            Tefferr = 'NULL'


# The log g
    puts('The ' + colored.yellow('logg'))
    logg = raw_input('> ')
    if logg == '':
        logg = 'NULL'
        loggerr = 'NULL'
    else:
        puts('The error on ' + colored.yellow('logg'))
        loggerr = raw_input('> ')
        if loggerr == '':
            loggerr = 'NULL'

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
    elif par.lower() == 'n' or par.lower() == 'no':
        puts(colored.magenta('\nCalculating the parallax...'))
        try:
            p = parallax(Teff, logg, V, M)
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
    vt = raw_input('> ')
    if vt == '':
        vt = 'NULL'
        vterr = 'NULL'
    else:
        puts('The error on '+colored.yellow('microturbulence'))
        vterr = raw_input('> ')
        if vterr == '':
            vt = 'NULL'

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

    zzz = ''
    with open(args.o, 'wb') as f:
        for z in params:
            zzz += str(z)
            zzz += '\t'
        zzz += 'NULL\n'
        f.write(zzz)
