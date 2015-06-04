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


def downloadExoplanet():
    """
    Download the table from exoplanetEU and save it to a file (exo.csv).

    Return a pandas DataFrame sorted in 'update'.
    """
    link_csv = 'http://www.exoplanet.eu/catalog/csv/?f=%22radial%22+IN'\
               'detection+OR+%22astrometry%22+IN+detection+OR+%22transit%22+'\
               'IN+detection'

    puts(colored.clean('Updating from exoplanets.eu.'))
    puts(colored.yellow('Please wait...'))
    response = urllib2.urlopen(link_csv)
    table = response.read()
    with open("exo.csv", "w") as f:
        f.write(table)

    exoplanet = pd.read_csv('exo.csv')
    return exoplanet


def sendingMail(names):
    """Sending the mail"""
    with open('mail.txt', 'r') as fp:
        msg = MIMEText(fp.read())

    with open('mailinfo.txt', 'r') as f:
        sender = f.readline().split(': ')[1].strip('\n')
        receiver = f.readline().split(': ')[1].strip('\n')
        smtp = f.readline().split(': ')[1].strip('\n')

    N = len(names)
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


if __name__ == '__main__':
    # Read the current version of SWEET-Cat
    names_ = ['name', 'hd', 'ra', 'dec', 'V', 'Verr', 'p', 'perr',
              'pflag', 'Teff', 'Tefferr', 'logg', 'logger',
              'n1', 'n2', 'vt', 'vterr', 'feh', 'feherr', 'M', 'Merr',
              'author', 'link', 'source', 'update', 'comment', 'n3']
    SC = pd.read_csv('WEBSITE_online.rdb', delimiter='\t', names=names_)
    sc_names = map(lambda x: x.lower().replace(' ', ''), SC.name)

    # Get all thene exoplanets from exoplanetEU (unique list)
    exoplanet = downloadExoplanet()
    exo_names = map(lambda x: x.lower().replace(' ', ''), exoplanet.star_name)

    NewStars = []
    for i, exo_name in enumerate(exo_names):
        if exo_name not in sc_names:
            new = exoplanet['# name'].values[i]
            NewStars.append(remove_planet(new))

    NewStars = list(set(NewStars))
    N = len(NewStars)

    if N:
        puts(colored.green(str(N) + " new exoplanet available!"))

        # Preparing list for SIMBAD
        simbad(NewStars, 'names.txt')

        # sendingMail(names)
    else:
        puts(colored.clean('No new updates seems to be available.'))
        puts(colored.clean('SWEET-Cat should be up to date'))
