import pandas as pd
import numpy as np
import warnings
from clint.textui import colored
warnings.simplefilter("ignore")


class Sweetcat:
    """Load SWEET-Cat database"""
    def __init__(self):

        # self.fname_sc = 'WEBSITE_online_EU-NASA_full_database.rdb'
        self.fname_sc = 'WEBSITE_online_EU-NASA_full_database_clean.rdb'

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

    # Select only the EU data
    sc_only_EU = sc.SC[sc.SC['database'].str.contains('EU')]

    # Drop the database column
    sc_EU = sc_only_EU.drop(columns=['database'])

    # Convert Tefferr column to integers
    sc_EU['Tefferr'] = sc_EU['Tefferr'].fillna('-111111')
    sc_EU['Tefferr'] = sc_EU['Tefferr'].astype(int).replace(-111111, 'NULL')

    # Replace NaN by NULL
    sc_EU.fillna(value='NULL', inplace=True)

    # Write to CSV file
    print('Write WEBSITE_online_EU.rdb')
    sc_EU.to_csv('WEBSITE_online_EU.rdb',
                 sep='\t', index=False, header=False)
