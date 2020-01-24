#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# My imports
import os
import urllib.request, urllib.error, urllib.parse
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
warnings.simplefilter("ignore")

def writeFile(fname, data):
    """Write data to a file"""
    with open(fname, 'w') as f:
        f.write(data)


class Update:
    """Check for updates to SWEET-Cat comparing with exoplanet.eu"""
    def __init__(self, controversial, download=False, nasa=False):
        self.controversial = controversial
        self.download = download
        self.nasa = nasa

        # if self.controversial:
        #     self.fname = 'exo_cont.csv'
        # else:
        if nasa:
            self.fname = 'nasaexo.csv'
        else:
            self.fname = 'exo.csv'

        self.blacklist = []
        # Kapteyn's can't be added with the ' in the website

        # Loading the SweetCat database
        self.readSC()

        # Loading the exoplanet database
        if self.nasa:
            self.downloadNasaExoplanet()
        else:
            self.downloadExoplanet()

    def downloadExoplanet(self):
        """
        Download the table from exoplanetEU and save it to a file (exo.csv).

        Return a pandas DataFrame sorted in 'update'.
        """
        if self.download:
            response = urllib.request.urlopen("http://exoplanet.eu/catalog/votable")
            table = response.read()
            with open('exo.xml', 'wb') as f:
                f.write(table)
            self.xml2csv()
        df = pd.read_csv(self.fname)
        df = df[(df.detection_type == 'Radial Velocity') | (df.detection_type == 'Primary Transit') | (df.detection_type == 'Astrometry')]
        self.exoplanet = df
        names=[self.remove_planet(x) for x in self.exoplanet['name']]
        self.exo_names = [x.strip() for x in names]

    def downloadNasaExoplanet(self):
        """
        Download the table from NASA exoplanet archive
        and save it to a file (nasaexo.csv).

        Return a pandas DataFrame sorted in 'update'.
        """

        urlRoot = "http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?"
        table = "&table=exoplanets"     # table of confirmed exoplanets
        select = "&select="             # which columns to select
        order = "&order=rowupdate"      # sort the table by date of last update
        tableFormat = "&format=csv"     # format of the table

        # Names of Sweet-Cat stellar columns
        starId = 'pl_hostname,hd_name,ra,dec,ra_str,dec_str,'
        detect = 'pl_name,pl_discmethod,'
        # keplerMag = 'st_optband,st_optmag,st_optmagerr,'
        vMag = 'st_vj,st_vjerr,'
        parallax = 'st_plx,st_plxerr1,st_plxerr2,'
        teff = 'st_teff,st_tefferr1,st_tefferr2,'
        logg = 'st_logg,st_loggerr1,st_loggerr2,'
        FeH = 'st_metfe,st_metfeerr1,st_metfeerr2,'
        mass = 'st_mass,st_masserr1,st_masserr2,'
        spt = 'st_spstr'

        # Choose SweetCat columns
        sweetCatColumns = starId + detect + vMag + parallax + teff + logg + FeH + mass + spt

        # Full url for the query
        nasa_url = urlRoot + table + select + sweetCatColumns + order + tableFormat

        # Download the data
        response = urllib.request.urlopen(nasa_url)
        table = response.read()

        # Write the NASA exoplanet archive
        with open('nasaexo.csv', 'wb') as f:
            f.write(table)

        # Load the NASA archive
        df = pd.read_csv(self.fname)
        # Select exoplanet detected with given methods
        df = df[(df.pl_discmethod == 'Radial Velocity') | (df.pl_discmethod == 'Transit') | (df.pl_discmethod == 'Astrometry')]
        self.exoplanet = df

        # List of names of the stars
        self.exo_names = [x.strip() for x in self.exoplanet['pl_hostname']]

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
            #os.remove('exo.xml')

    def remove_planet(self, name):
        """Remove the trailing b, c, d, etc in the stellar name"""
        planets = 'bcdefghijB'
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
                  'author', 'link', 'source', 'update', 'comment', 'database',
                  'n3']
        
        # SC = pd.read_csv('WEBSITE_online.rdb', delimiter='\t', names=names_)
        SC = pd.read_csv('website_nasa-eu_updated.rdb',
                         delimiter='\t', names=names_)

        self.sc_names = [x.lower().replace(' ', '').replace('-', '') for x in SC.name]
        self.sc_names = list(map(str.strip, self.sc_names))
        self.sc_names_orig = [x.strip() for x in SC.name]
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
        """
        Check the difference the exoplanet database (EU or NASA)
        and the SweetCat database
        INPUTS: self = exoplanet database in pandas dataframe
        OUTPUTS: names.txt = file with stars that are not in SweetCat but
                             stars are in the exoplanet database (EU or NASA)
        Prints the stars that are in SweetCat
        and that are not in the exoplanet database (EU or NASA)
        """

        # We have this already, but without the ' in the name.
        print('\n    Matching data base...')

        NewStars = []
        coordExo = coord.SkyCoord(ra=self.exoplanet['ra'].values,
                                  dec=self.exoplanet['dec'].values,
                                  unit=(u.deg, u.deg),
                                  frame='icrs')
        coordSC = coord.SkyCoord(ra=self.coordinates['ra'].values,
                                 dec=self.coordinates['dec'].values,
                                 unit=(u.hourangle, u.deg),
                                 frame='icrs')

        for i, exo_name in enumerate(self.exo_names):
            starName = exo_name
            # starName = self.exoplanet['star_name'].values[i]

            # Clean star name
            tmp = starName.lower().replace(' ', '').replace('-', '')

            sep = coordExo[i].separation(coordSC).arcsecond
            ind = np.where(sep < 5.)[0]

            if len(ind) == 0:
                try:
                    # it didn't find by position but it finds by name
                    position = self.sc_names.index(tmp)

                except:
                    position = -1
                    # it didn't find by position and neither by name
                    if (tmp not in self.blacklist):
                        NewStars.append(starName)


                #  REMOVE THE COMMENTS TO CHECK PROBLEMS IN POSITION
#                 if position>=0:

#                     result_table = Simbad.query_object(new)
#                     try:
#                         # check in Simbad the position and see where the coordinates are wrong
#                         ra=str(result_table['RA'][0])
#                         dec=str(result_table['DEC'][0])
#                         coordS=coord.SkyCoord(ra, dec, unit=(u.hourangle, u.degree), frame='icrs')
#                         sepES = coordExo[i].separation(coordS).arcsecond
#                         sepSCS = coordSC[position].separation(coordS).arcsecond

#                         if sepES>sepSCS:
# #                            print new,': has wrong position in Exoplanet.eu, it is: ',coordExo[i].ra.deg, coordExo[i].dec.deg,', but should be: ',coordS.ra.deg,coordS.dec.deg
#                             pass
#                         else:
#                            print new,': has wrong position in Sweet-Cat, it is: ',coordSC[position].ra.deg, coordSC[position].dec.deg,', but should be: ',coordS.ra.deg,coordS.dec.deg
#                     except:
#                         print 'Star not found in Simbad with this name ',new,'.\n Position in Exoplanet.eu:',coordExo[i].ra.deg, coordExo[i].dec.deg,'\n Position in SC: ',coordSC[position].ra.deg, coordSC[position].dec.deg
#                         with open('starnotfoundinsimbad.list', 'a') as f:
#                             f.write(new+'\n')

            #  REMOVE THE COMMENTS TO CHECK PROBLEMS WITH THE NAMES
            # else:
            #     position=ind[0]
            #     if (tmp not in self.sc_names):
            #         print new, ': has a different name in Sweet-Cat (',self.sc_names[position],')'

        NewStars = sorted(list(set(NewStars)))
        Nstars = len(NewStars)

        if Nstars:
            puts('    '+colored.green(str(Nstars) + " new exoplanet available!"))
            writeFile('names.txt', '\n'.join(NewStars))
            updated=False
        else:
            puts(colored.clean('    No new updates available.'))
            updated=True

        # removing planets that are not in Exoplanet.eu anymore    
        NewStars = []

        for i, scname in enumerate(self.sc_names_orig):
            sep = coordSC[i].separation(coordExo).arcsecond
            ind=np.where(sep<5.)[0]
            if len(ind)==0:
                try:
                    # it didn't find by position but it finds by name
                    position=self.exo_names.index(scname)
                except:
                    position=-1
                    # it didn't find by position and neither by name
                    if (tmp not in self.blacklist):
                        NewStars.append(scname)

        NewStars = sorted(list(set(NewStars)))
        Nstars = len(NewStars)
        if Nstars:
            puts(colored.green('    '+str(Nstars) + " exoplanet has to be removed!"))
            print('\n    '.join(NewStars))
        else:
            puts(colored.clean('    No planet to remove.'))
            if updated:        
                puts(colored.clean('    SWEET-Cat is up to date'))
                puts(colored.green('    Great job :)'))


if __name__ == '__main__':
    with open('starnotfoundinsimbad.list', 'a') as f:
        f.write(str(time.strftime("%d-%m-%Y"))+'\n')
    exo_database = Update(controversial=False, download=True, nasa=False)
    exo_database.update()
