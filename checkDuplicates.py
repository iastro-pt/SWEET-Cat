import pandas as pd
from astropy import coordinates as coord
from astropy import units as u


class Sweetcat:
    """Load SWEET-Cat database"""
    def __init__(self):

        self.fname_sc = 'WEBSITE_online_EU-NASA_full_database.rdb'
        
        # Loading the SweetCat database
        self.readSC()

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


if __name__ == '__main__':

    # Loading SWEET Cat
    sc = Sweetcat()

    # Check for duplicates, subset of columns can be changed
    print(sc.SC[sc.SC.duplicated(['ra', 'dec'], keep=False)])

    # Indexes of the duplicates
    indexes = sc.SC[sc.SC.duplicated(['ra', 'dec'], keep=False)].index

    # Remove a row
    # new_sc = sc.SC.drop([2728])
    # new_sc.to_csv('WEBSITE_online_EU-NASA_full_database_minusHD21749.rdb',
    #               sep='\t', index=False, header=False)

    # Select only the EU data
    sc_EU = new_sc[new_sc['database'].str.contains('EU')]
    # Drop the database column
    sc_like_old = sc_EU.drop(columns=['database'])
    #sc_like_old.to_csv('WEBSITE_online_EU-updated_04-03-2020.rdb',
    #                   sep='\t', index=False, header=False)
