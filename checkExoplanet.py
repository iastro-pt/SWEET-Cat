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

        self.blacklist = ['HD 59686 A','Kepler-420 A','Kepler-539','HD 202206 B']

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
        planets = 'abcdefgh'
        for planet in planets:
            if name.endswith(' %s' % planet):
                return name[:-2]
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

    def _sccoordinates(self, idx):
        '''Transform the coordinates to degrees

        Input
        -----
        idx : int
          The index of the SWEET-Cat catalogue

        Output
        ------
        RAsc : float
          RA in degrees
        DEsc : float
          DEC in degrees
        '''
        aux = self.coordinates['ra'].values[idx]
        RAsc = (float(aux[0:2])+float(aux[3:5])/60.+float(aux[6:])/3600.)*15.

        aux = self.coordinates['dec'].values[idx]
        if aux[0] == '-':
            DEsc = float(aux[0:3])-float(aux[4:6])/60.-float(aux[7:])/3600.
        else:
            DEsc = float(aux[0:3])+float(aux[4:6])/60.+float(aux[7:])/3600.
        return RAsc, DEsc

    def update(self):
        if self.controversial:
            df = pd.read_csv('exo.csv')
            true_names = map(lambda x: self.remove_planet(x).lower().replace(' ', ''), df.name)

        # We have this already, but without the ' in the name.
        print '\nMatching data base...'
        NewStars = []
        for i, exo_name in enumerate(self.exo_names):
            new = self.remove_planet(self.exoplanet['name'].values[i])
            tmp = new.lower().replace(' ', '')

            RAexo = self.exoplanet['ra'].values[i]
            DEexo = self.exoplanet['dec'].values[i]
            found = False

            for j in xrange(len(self.sc_names)):
                # converting RA (h:m:s->degrees) and DE (d:m:s->degrees)
                RAsc, DEsc = self._sccoordinates(j)

                # compute the spherical distance
                dist = 999
                if abs(RAexo-RAsc) < 1:
                    dist = 3600 * (((RAexo-RAsc)*np.cos(np.radians(DEexo)))**2 + (DEexo-DEsc)**2)**0.5
                if dist < 5:  # 5 arc seconds is a good tolerance in distance
                    found = True
                    break

            if self.controversial:
                if exo_name in self.sc_names:
                    if (tmp in self.sc_names) and (tmp not in true_names):
                        NewStars.append(new)
            else:
                if (exo_name not in self.sc_names) and (exo_name not in self.blacklist):
                    if new in self.blacklist:
                        continue

                    # the positions in exoplanets.eu are wrong sometimes because they dont use J2000, so double check position and name
                    elif (not found) and (tmp not in self.sc_names):
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
