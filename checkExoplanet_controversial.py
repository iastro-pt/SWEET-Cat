#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import print_function
import urllib2
from datetime import datetime
import pandas as pd
import numpy as np

from Simbad import simbad

# For fun, but still useful
import smtplib
from email.mime.text import MIMEText
from clint.textui import puts, colored


def planetInString(planets):
    """
    Simbad does not like planet names, so remove 'b', 'c' etc. in the of the
    star name
    """
    for i, planet in enumerate(planets):
        for letter in 'bcdefg':
            if planet.endswith(' %s' % letter):
                planets[i] = planet[0:-2]
    return planets


def downloadExoplanet():
    link_csv = 'http://exoplanet.eu/catalog/csv/?f=%22controversial'\
               '%22+IN+detection'

    puts(colored.clean('Updating from exoplanets.eu.'))
    puts(colored.yellow('Please wait...'))
    response = urllib2.urlopen(link_csv)
    table = response.read()
    with open("exo_cont.csv", "w") as f:
        f.write(table)

    exoplanet = pd.read_csv('exo.csv')
    exoplanet.updated = pd.to_datetime(exoplanet.updated)
    exoplanet.sort(['updated'], ascending=False, inplace=True)
    return exoplanet


def update(exoplanet, starsID, update_sweetcat):
    return exoplanet['# name'][exoplanet.updated > update_sweetcat]


def sendingMail(names):
    """Sending the mail"""
    with open('mail.txt', 'r') as fp:
        msg = MIMEText(fp.read())

    with open('mailinfo.txt', 'r') as f:
        sender = f.readline().split(': ')[1].strip('\n')
        receiver = f.readline().split(': ')[1].strip('\n')
        smtp = f.readline().split(': ')[1].strip('\n')

    msg['Subject'] = 'Update available to SWEET-Cat: %s new controversial'\
                     ' exoplanets' % N
    msg['From'] = sender
    msg['To'] = receiver

    try:
        s = smtplib.SMTP(smtp)
        s.sendmail(sender, [receiver], msg.as_string())
        s.quit()
    except smtplib.SMTPRecipientsRefused:
        puts(colored.red('Not able to send an email...'))
        raise smtplib.SMTPRecipientsRefused


if __name__ == '__main__':
    # Read the current version of SWEET-Cat
    names_ = ['name', 'hd', 'ra', 'dec', 'V', 'Verr', 'p', 'perr',
              'pflag', 'Teff', 'Tefferr', 'logg', 'logger',
              'n1', 'n2', 'vt', 'vterr', 'feh', 'feherr', 'M', 'Merr',
              'author', 'link', 'source', 'update', 'comment', 'n3']
    SC = pd.read_csv('WEBSITE_online.rdb', delimiter='\t', names=names_)
    starsID = SC.name  # All the stars from SWEET-Cat

    # Find the latest update time
    newest = max(SC.values[:, 24]).replace('-', ' ')
    update_sweetcat = datetime.strptime(newest, '%Y %m %d')

    # Get all the new exoplanet after the newest update time
    exoplanet = downloadExoplanet()
    names = update(exoplanet, starsID, update_sweetcat)
    names = np.array(list(set(planetInString(names.values))))
    # Remove updated exoplanets, and only look at unique hosts
    names = [name for name in names if not len(SC[SC.name == name].values)]
    N = len(names)

    if N:
        puts(colored.green(str(N) + " new exoplanet available!"))

        # Preparing list for SIMBAD
        simbad(names, 'names_controversial.txt')

        sendingMail(names)
    else:
        puts(colored.clean('No new updates seems to be available.'))
        puts(colored.clean('SWEET-Cat should be up to date'))
