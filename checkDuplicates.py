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

    # Check for duplicates, subset of columns can be changed
    print('\nChecking for possible duplicates ...')
    print(colored.green('Same RA/DEC'))
    print(sc.SC[sc.SC.duplicated(['ra', 'dec'], keep=False)][['name',
                                                              'hd',
                                                              'ra',
                                                              'dec']])
    print(colored.green('\nSame HD number'))
    print(sc.SC[sc.SC.duplicated(['hd'],
                                 keep=False)].dropna(subset=['hd'])[['name',
                                                                     'hd',
                                                                     'ra',
                                                                     'dec']])
    print(colored.green('\nApproximate RA/DEC ...'))
    # Remove the characters after the . in the coordinates
    ra_sc = sc.SC['ra'].values.tolist()
    ra_approx = list(map(lambda i: i[:i.find('.')], ra_sc))
    dec_sc = sc.SC['dec'].values.tolist()
    dec_approx = list(map(lambda i: i[:i.find('.')], dec_sc))

    # Check for similar RA/DEC
    idx_duplicate = []
    for idx, (ra, dec) in enumerate(zip(ra_approx, dec_approx)):
        dupli = list(np.where((np.array(ra_approx) == ra) &
                              (np.array(dec_approx) == dec))[0])
        if len(dupli) > 1:
            idx_duplicate.append(dupli)

    # Print possible duplicates
    print(colored.green('RA/DEC are similar: possible duplicates\n'))
    unique_duplicate = set([tuple(t) for t in idx_duplicate])
    for idx in unique_duplicate:
        print(sc.SC.iloc[list(idx)][['name', 'hd', 'ra', 'dec']])


    # Uncomment some of the following lines to remove duplicates
    # Indexes of the duplicates
    # indexes = sc.SC[sc.SC.duplicated(['ra', 'dec'], keep=False)].index

    # Remove a row (HD21749)
    # new_sc = sc.SC.drop([2728])


    # # Write the new file
    # # Convert Tefferr column to integers
    # new_sc['Tefferr'] = new_sc['Tefferr'].fillna('-111111')
    # new_sc['Tefferr'] = new_sc['Tefferr'].astype(int).replace(-111111, 'NULL')

    # # Replace NaN by NULL
    # new_sc.fillna(value='NULL', inplace=True)
    # new_sc.to_csv('WEBSITE_online_EU-NASA_full_database_clean_minusKELT-22A.rdb',
    #               sep='\t', index=False, header=False)

    # Select only the EU data
    # sc_EU = new_sc[new_sc['database'].str.contains('EU')]
    # Drop the database column
    # sc_like_old = sc_EU.drop(columns=['database'])
    #sc_like_old.to_csv('WEBSITE_online_EU-updated_04-03-2020.rdb',
    #                   sep='\t', index=False, header=False)
