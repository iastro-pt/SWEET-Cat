#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
import os
import urllib2
import pandas as pd
import numpy as np
from astropy.io import votable
import warnings

# For fun, but still useful
from clint.textui import puts, colored


def writeFile(fname, data):
    """Write data to a file"""
    with open(fname, 'w') as f:
        f.write(data)


class Update:
    """Check for updates to SWEET-Cat comparing with exoplanet.eu"""
    def __init__(self, controversial, download=False):
        self.controversial = controversial
        self.download = download

        if self.controversial:
            self.fname = 'exo_cont.csv'
        else:
            self.fname = 'exo.csv'

        self.blacklist = ['Kapteyn\'s', 'KELT-14', 'K2-8', 'K2-16', 'K2-19',
                          'KOI-2939', 'KOI-2828', 'Kepler-102', 'KIC-10024862',
                          'KIC-8012732', 'KIC-9413313', 'EPIC 388', 'BD+14 4559',
                          'BD+20 2457', 'HD 47 32', 'GJ 86 A', 'HD 4732',
                          'KOI-2939 (AB)', 'Kepler-453 (AB)', 'Kepler-64 (AB)',
                          'Kepler-539', 'Kepler-1647 (AB)', 'K2-24', 'CVSO 30 b ',
                          'HD 59686 A', 'HAT-47', 'HAT-P-27-WASP-40',
                          'HAT-P-30-WASP-51','GJ 221  BD-06 1339','K2-99','K2-33','KOI-1089.02','KOI-1299','KOI-368.01','KOI-4427.01','Kepler-13 A']

#,'eta Cnc','eta Umi','ups An','55 Cn' ***star with letters missing***

        self.readSC()
        self.downloadExoplanet()

    def downloadExoplanet(self):
        """
        Download the table from exoplanetEU and save it to a file (exo.csv).

        Return a pandas DataFrame sorted in 'update'.
        """
        if self.download:
            response = urllib2.urlopen("http://exoplanet.eu/catalog/votable")
            table = response.read()
            with open('exo.xml', 'w') as f:
                f.write(table)
            self.xml2csv()
        df = pd.read_csv(self.fname)
        df = df[(df.detection_type == 'Radial Velocity') | (df.detection_type == 'Primary Transit') | (df.detection_type == 'Astrometry')]
        self.exoplanet = df
        self.exo_names = self.exoplanet.star_name

    def xml2csv(self):
        """Convert the saved xml file to csv and read with pandas"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            vo = votable.parse('exo.xml', invalid='mask', pedantic=False)
            vo = vo.get_first_table().to_table(use_names_over_ids=True)
            df = vo.to_pandas()

            # Divide the data in Confirmed and not.
            df[df.planet_status == 'Confirmed'].to_csv('exo.csv', index=False)
            df[(df.planet_status == 'Unconfirmed') | (df.planet_status == 'Candidate')].to_csv('exo_cont.csv', index=False)
            os.remove('exo.xml')

    def remove_planet(self, name):
        """Remove the trailing b, c, d, etc in the stellar name"""
        planets = ' a, b, c, d, e, f, g, h'.split(',')
        for planet in planets:  # Probably not more planets currently
            if name.endswith(planet):
                return name.strip(' %s' % planet)
        return name

    def readSC(self):
        # TODO: Use the ra and dec, and match with coordinates instead of name
        #       stored in self.coordinates.

        # Read the current version of SWEET-Cat
        names_ = ['name', 'hd', 'ra', 'dec', 'V', 'Verr', 'p', 'perr',
                  'pflag', 'Teff', 'Tefferr', 'logg', 'logger',
                  'n1', 'n2', 'vt', 'vterr', 'feh', 'feherr', 'M', 'Merr',
                  'author', 'link', 'source', 'update', 'comment', 'n3']
        SC = pd.read_csv('WEBSITE_online.rdb', delimiter='\t', names=names_)
        self.sc_names = map(lambda x: x.lower().replace(' ', ''), SC.name)
        self.sc_names = map(str.strip, self.sc_names)
        self.coordinates = SC.loc[:, ['ra', 'dec']]

    def update(self):
        if self.controversial:
            df = pd.read_csv('exo.csv')
            true_names = map(lambda x: self.remove_planet(x.lower().replace(' ', '')), df.name)

        # We have this already, but without the ' in the name.
        NewStars = []
        for i, exo_name in enumerate(self.exo_names):
            new = self.remove_planet(self.exoplanet['name'].values[i])
            tmp = new.lower().replace(' ', '')

            if self.controversial:
                if exo_name in self.sc_names:
                    if (tmp in self.sc_names) and (tmp not in true_names):
                        NewStars.append(new)
            else:
                if (exo_name not in self.sc_names) and (exo_name not in self.blacklist):
                    if new in self.blacklist:
                        continue
                    if tmp not in self.sc_names:
                        NewStars.append(new)

        NewStars = sorted(list(set(NewStars)))
        Nstars = len(NewStars)

        if Nstars:
            puts(colored.green(str(Nstars) + " new exoplanet available!"))

            if self.controversial:
                writeFile('names_contr.txt', '\n'.join(NewStars))
            else:
                writeFile('names.txt', '\n'.join(NewStars))
        else:
            puts(colored.clean('No new updates available.'))
            puts(colored.clean('SWEET-Cat is up to date'))
            puts(colored.green('Great job :)'))


if __name__ == '__main__':
    new = Update(controversial=False, download=True)
    new.update()
    print('\n')
    remove = Update(controversial=True)
    remove.update()
