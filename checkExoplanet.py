#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
import os
import urllib2
import pandas as pd
import numpy as np
import time
from astropy.io import votable
import warnings
from astropy import coordinates as coord
from astropy import units as u
from astroquery.simbad import Simbad

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

        self.blacklist = ['Kepler-1647 (AB)','Kepler-453 (AB)']
# Kepler-1647 (AB), Kepler-453 (AB) keep in the black list because I have the parameters for both stars
# Kapteyn's can't be added with the ' in the website

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
        planets = 'bcdefghB'
        for planet in planets:
            if name.endswith(' %s' % planet):
                return name[:-2]
        # some exoplanets have .01 or .02 in the name 
        if name.endswith('.01') or name.endswith('.02') or name.endswith('.2'):
            return name[:-3]                
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
        self.sc_names = map(lambda x: x.lower().replace(' ', '').replace('-', ''), SC.name)
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
            starName=self.exoplanet['star_name'].values[i]
            tmp = new.lower().replace(' ', '').replace('-', '') 

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
                    position=j
                    break
            if self.controversial:
                if tmp in self.sc_names and (tmp not in true_names):
                    NewStars.append(new)
            else:
                if (not found):
                    try:
                        # it didn't find by position but it finds by name
                        position=self.sc_names.index(tmp)
                    except:
                        position=-1
                        # it didn't find by position and neither by name
                        if (tmp not in self.blacklist):
                            NewStars.append(new)

                    #  REMOVE THE COMMENTS TO CHECK PROBLEMS IN POSITION
#                    if position>=0:
#                        RAsc, DEsc = self._sccoordinates(position)
#                        result_table = Simbad.query_object(new)
                        # it found the star in Simbad
#                        try:
                            # check in Simbad the position and see where the coordinates are wrong
#                            ra=result_table['RA'][0]
#                            dec=result_table['DEC'][0]
                            # position in Simbad
#                            coordS=coord.SkyCoord(ra, dec, unit=(u.hourangle, u.degree), frame='icrs')
                            # position in Exoplanet.eu 
#                            coordE = coord.SkyCoord(RAexo, DEexo, unit=(u.degree, u.degree), frame='icrs')
                            # position in Sweet-Cat
#                            coordSC = coord.SkyCoord(RAsc, DEsc, unit=(u.degree, u.degree), frame='icrs')
#                            sepES = coordE.separation(coordS).arcsecond
#                            sepSCS = coordS.separation(coordSC).arcsecond
#                            if sepES>sepSCS:
#                                print new,': has wrong position in Exoplanet.eu, it is: ',coordE.ra, coordE.dec,', but should be: ',coordS.ra,coordS.dec
#                            else:
#                                print new,': has wrong position in Sweet-Cat, it is: ',coordSC.ra, coordSC.dec,', but should be: ',coordS.ra,coordS.dec
#                        except:    
#                            with open('starnotfoundinsimbad.list', 'a') as f:
#                                f.write(new+'\n')

                #  REMOVE THE COMMENTS TO CHECK PROBLEMS WITH THE NAMES
#                else:
#                    if (tmp not in self.sc_names):
#                        print new, ': has a different name in Sweet-Cat (',self.sc_names[position],')'
                        
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
    with open('starnotfoundinsimbad.list', 'a') as f:
        f.write(str(time.strftime("%d-%m-%Y"))+'\n')
    new = Update(controversial=False, download=True)
    new.update()
    print('\n')
    remove = Update(controversial=True)
    remove.update()
