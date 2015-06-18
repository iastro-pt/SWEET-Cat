#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
# from __future__ import print_function
import urllib2
from datetime import datetime
import pandas as pd
import numpy as np
from Simbad import simbad

# For fun, but still useful
import smtplib
from email.mime.text import MIMEText
from clint.textui import puts, colored


def downloadExoplanet(link, controversial=False):
    """
    Download the table from exoplanetEU and save it to a file (exo.csv).

    Return a pandas DataFrame sorted in 'update'.
    """

    puts(colored.clean('Updating from exoplanets.eu.'))
    puts(colored.yellow('Please wait...'))
    response = urllib2.urlopen(link)
    table = response.read()
    if controversial:
        with open("exo_cont.csv", "w") as f:
            f.write(table)
        return pd.read_csv('exo_cont.csv')
    else:
        with open("exo.csv", "w") as f:
            f.write(table)

        return pd.read_csv('exo.csv')


def sendingMail(names, controversial=False):
    """Sending the mail"""
    with open('mail.txt', 'r') as fp:
        msg = MIMEText(fp.read())

    with open('mailinfo.txt', 'r') as f:
        sender = f.readline().split(': ')[1].strip('\n')
        receiver = f.readline().split(': ')[1].strip('\n')
        smtp = f.readline().split(': ')[1].strip('\n')

    N = len(names)
    if controversial:
        msg['Subject'] = 'Update available to SWEET-Cat: %i new controversial'\
                         ' exoplanets' % N
    else:
        msg['Subject'] = 'Update available to SWEET-Cat: %i new exoplanets' % N
    msg['From'] = sender
    msg['To'] = receiver

    try:
        s = smtplib.SMTP(smtp)
        s.sendmail(sender, [receiver], msg.as_string())
        s.quit()
    except smtplib.SMTPRecipientsRefused:
        puts(colored.red('Not able to send an email...'))
        raise smtplib.SMTPRecipientsRefused


def remove_planet(name):
    """Remove the trailing b, c, d, etc in the stellar name"""
    for planet in 'abcdefgh':  # Probably not more planets currently
        if name.endswith(planet):
            return name.strip(' %s' % planet)
    return name


def main(link, controversial=False):
    # Read the current version of SWEET-Cat
    names_ = ['name', 'hd', 'ra', 'dec', 'V', 'Verr', 'p', 'perr',
              'pflag', 'Teff', 'Tefferr', 'logg', 'logger',
              'n1', 'n2', 'vt', 'vterr', 'feh', 'feherr', 'M', 'Merr',
              'author', 'link', 'source', 'update', 'comment', 'n3']
    SC = pd.read_csv('WEBSITE_online.rdb', delimiter='\t', names=names_)
    sc_names = map(lambda x: x.lower().replace(' ', ''), SC.name)
    sc_names = map(str.strip, sc_names)

    # Get all thene exoplanets from exoplanetEU (unique list)
    exoplanet = downloadExoplanet(link, controversial=controversial)
    exo_names = map(lambda x: x.lower().replace(' ', ''), exoplanet.star_name)
    if controversial:
        df = pd.read_csv('exo.csv')
        true_names = map(lambda x: remove_planet(x.lower().replace(' ', '')), df.star_name)

    # We have this already, but without the ' in the name.
    blacklist = ['Kapteyn\'s']
    NewStars = []
    for i, exo_name in enumerate(exo_names):
        if controversial:
            if exo_name in sc_names:
                new = remove_planet(exoplanet['# name'].values[i])
                tmp = new.lower().replace(' ', '')
                if tmp in sc_names and tmp not in true_names:
                    NewStars.append(new)
        else:
            if exo_name not in sc_names:
                new = remove_planet(exoplanet['# name'].values[i])
                if new in blacklist:
                    continue
                if new.lower().replace(' ', '') not in sc_names:
                    NewStars.append(new)

    NewStars = sorted(list(set(NewStars)))
    N = len(NewStars)

    if N:
        puts(colored.green(str(N) + " new exoplanet available!"))

        # Preparing list for SIMBAD
        if controversial:
            simbad(NewStars, 'names_contr.txt')
        else:
            simbad(NewStars, 'names.txt')

        # sendingMail(names, controversial=controversial)
    else:
        puts(colored.clean('No new updates seems to be available.'))
        puts(colored.clean('SWEET-Cat should be up to date'))


if __name__ == '__main__':
    link1 = 'http://www.exoplanet.eu/catalog/csv/?f=%22radial%22+IN'\
            'detection+OR+%22astrometry%22+IN+detection+OR+%22transit%22+'\
            'IN+detection'
    link2 = 'http://exoplanet.eu/catalog/csv/?f=%22controversial'\
            '%22+IN+detection'

    # main(link1, controversial=False)
    # print('\n')
    main(link2, controversial=True)
