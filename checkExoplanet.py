#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# My imports
import os, sys
import time
import argparse
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


def _parse_args():
    parser = argparse.ArgumentParser('Check databases for new exoplanets')
    parser.add_argument('database', choices=['EU', 'NASA'],
                        type=str, help='Which database to check (EU or NASA)')
    parser.add_argument('-d', '--download', action='store_true',
                        help='Download the latest database')
    args = parser.parse_args()
    return args

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
        dir_script = os.path.dirname(__file__)
        if dir_script != "" :
            dir_script+='/'
        self.fname_sc = dir_script+'WEBSITE_online_EU-NASA_full_database_clean.rdb'

        # if self.controversial:
        #     self.fname = 'exo_cont.csv'
        # else:
        if nasa:
            self.fname = dir_script+'nasaexo.csv'
        else:
            self.fname = dir_script+'exo.csv'

        self.blacklist = ['hatp69', 'hatp70', 'denispj082303.1491201']
        # HAT stars are A stars and denis is a L dwarf
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
        names = [self.remove_planet(x) for x in self.exoplanet['name']]
        self.exo_names = [x.strip() for x in names]
        self.exo_names_clean = [x.lower().replace(' ', '').replace('-', '')
                                for x in names]
        self.exo_names_clean = list(map(str.strip, self.exo_names_clean))

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
        self.exo_names_clean = [x.lower().replace(' ', '').replace('-', '')
                                for x in self.exoplanet['pl_hostname']]
        self.exo_names_clean = list(map(str.strip,
                                        self.exo_names_clean))

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
        SC = pd.read_csv(self.fname_sc, delimiter='\t', names=names_)

        # Clean star names
        self.sc_names = [x.lower().replace(' ', '').replace('-', '') for x in SC.name]
        self.sc_names = list(map(str.strip, self.sc_names))
        # Original star names
        self.sc_names_orig = [x.strip() for x in SC.name]
        # Coordinates of the stars in SWEET-Cat
        self.coordinates = SC.loc[:, ['ra', 'dec']]
        # SWEET-Cat (used to automatically update the database label)
        self.SC = SC

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
                 fname_sc_day-month_hour:minute = new SWEET-Cat dabatase
                                                  with updated database label

        Prints the stars that are in SweetCat
        and that are not in the exoplanet database (EU or NASA)
        """

        # We have this already, but without the ' in the name.
        print('\n    Matching database ...')
        if self.nasa:
            print('    NASA exoplanet archive')
        else:
            print('    Extrasolar Planets Encyclopaedia')

        NewStars = []
        # Star corrdinates in the exoplanet database
        coordExo = coord.SkyCoord(ra=self.exoplanet['ra'].values,
                                  dec=self.exoplanet['dec'].values,
                                  unit=(u.deg, u.deg),
                                  frame='icrs')
        # Star coordinates in the SWEET-Cat database
        coordSC = coord.SkyCoord(ra=self.coordinates['ra'].values,
                                 dec=self.coordinates['dec'].values,
                                 unit=(u.hourangle, u.deg),
                                 frame='icrs')

        # -------------------------------------------------------
        # Adding stars that are not
        # in Exoplanet.EU nor NASA exoplanet archive
        # For stars which already are in one of the database,
        # only the database label is updated by this script
        # -------------------------------------------------------
        for i, exo_name in enumerate(self.exo_names):
            starName = exo_name
            # starIndex = self.exoplanet[self.exoplanet['pl_hostname'] == exo_name].index[0]
            # starName = self.exoplanet['star_name'].values[i]

            # Clean star name
            tmp = starName.lower().replace(' ', '').replace('-', '')

            # Check if the star in exoplanet database is in SWEET-Cat
            sep = coordExo[i].separation(coordSC).arcsecond
            ind = np.where(sep < 5.)[0]

            # Star is already in SWEET-Cat
            if len(ind) != 0.:
                ind_SC = ind[0]

                # Check if the name of the database (EU and/or NASA)
                # is written in the SWEET-CAT file, if not it's added
                if self.nasa:
                    if 'NASA' in self.SC.loc[ind_SC].database:
                        pass
                    else:
                        print('\nChecking star: ', starName, '(found by position)')
                        # print(self.exoplanet.loc[starIndex][['pl_hostname', 'ra_str', 'dec_str']])
                        # print(self.SC.loc[ind_SC][['name', 'ra', 'dec', 'database', 'n3']])
                        print(' > adding NASA label')
                        self.SC.at[ind_SC, 'database'] = self.SC.at[ind_SC, 'database'] + ',NASA'
                else:
                    if 'EU' in self.SC.loc[ind_SC].database:
                        pass
                    else:
                        print('\nChecking star: ', starName, '(found by position)')
                        print(' > adding EU label')
                        self.SC.at[ind_SC, 'database'] = self.SC.at[ind_SC, 'database'] + ',EU'

            if len(ind) == 0:
                try:
                    # it didn't find by position but it finds by name
                    position = self.sc_names.index(tmp)

                    # Check the name of the database (EU and/or NASA)
                    if self.nasa:
                        if 'NASA' in self.SC.loc[position].database:
                            pass
                        else:
                            print('\nChecking star: ', starName, i, '(found by name)')
                            print(' > adding NASA label')
                            self.SC.at[position, 'database'] = self.SC.at[position, 'database'] + ',NASA'
                    else:
                        if 'EU' in self.SC.loc[position].database:
                            pass
                        else:
                            print('\nChecking star: ', starName, i, '(found by name)')
                            print(' > adding EU label')
                            self.SC.at[position, 'database'] = self.SC.at[position, 'database'] + ',EU'

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
            updated = False
        else:
            puts(colored.clean('\n    No new updates available.'))
            updated = True

        # -------------------------------------------------------
        # Removing planets that are not
        # in Exoplanet.EU and NASA exoplanet archive anymore
        # -------------------------------------------------------
        print('\n    Checking for stars to remove ...')
        NewStars = []

        for i, scname in enumerate(self.sc_names_orig):
            # Clean star name
            tmp = starName.lower().replace(' ', '').replace('-', '')

            # Check if the star in SWEET-Cat is in exoplanet database
            sep = coordSC[i].separation(coordExo).arcsecond
            ind = np.where(sep < 5.)[0]

            if len(ind) == 0:
                try:
                    # it didn't find by position but it finds by name
                    position = self.exo_names_clean.index(tmp)

                except:
                    # Star in Sweet-Cat is not found
                    # in the exoplanet database (EU or NASA)
                    position = -1

                    # Check if the star is not from the other database
                    if self.nasa:
                        if 'EU' in self.SC.loc[i].database:
                            # Star is in EU
                            # Star will not be removed from SWEET-Cat
                            if 'NASA' in self.SC.loc[i].database and scname != 'Barnards':
                                # Removing NASA label from database column
                                print(' > Removing NASA label')
                                self.SC.at[i, 'database'] = self.SC.at[i, 'database'].replace('NASA', '').replace(',', '')
                            continue
                    else:
                        if 'NASA' in self.SC.loc[i].database:
                            # Star is in NASA
                            # Star will not be removed from SWEET-Cat
                            if 'EU' in self.SC.loc[i].database and scname != 'Barnards':
                                # Removing EU label from database column
                                print(scname)
                                print(' > Removing EU label')
                                self.SC.at[i, 'database'] = self.SC.at[i, 'database'].replace('EU', '').replace(',', '')
                            continue

                    # it didn't find by position and neither by name
                    # star is not from the other database
                    if (tmp not in self.blacklist):
                        NewStars.append(scname)

        NewStars = sorted(list(set(NewStars)))
        Nstars = len(NewStars)
        if Nstars:
            puts(colored.green('    '+str(Nstars) + " exoplanet has to be removed!"))
            print('\n    '.join(NewStars))
        else:
            puts(colored.clean('\n    No star to remove.'))
            if updated:
                puts(colored.clean('    SWEET-Cat is up to date'))
                puts(colored.green('    Great job :)'))

        # Date and time
        timestr = time.strftime("%d-%m-%H:%M")
        filename = os.path.splitext(self.fname_sc)[0] + '_' + timestr + '.rdb'

        # Convert Tefferr column to integers
        self.SC['Tefferr'] = self.SC['Tefferr'].fillna('-111111')
        self.SC['Tefferr'] = self.SC['Tefferr'].astype(int).replace(-111111,
                                                                    'NULL')
        # Replace NaN by NULL
        self.SC.fillna(value='NULL', inplace=True)

        # Write new SWEET-Cat database
        print('\n    Writing the file: ', filename)
        self.SC.to_csv(filename, sep='\t', index=False, header=False)


if __name__ == '__main__':
    args = _parse_args()
    print(args)

    with open('starnotfoundinsimbad.list', 'a') as f:
        f.write(str(time.strftime("%d-%m-%Y"))+'\n')

    # Load SWEET-CAT and EU/NASA databases
    nasa = args.database == 'NASA'
    exo_database = Update(controversial=False, download=args.download, nasa=nasa)
    # Check for new planet host stars and add the names of EU/NASA databases
    exo_database.update()
