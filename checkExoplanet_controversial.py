#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
import csv
import operator
import urllib2
from datetime import datetime

from Simbad import simbad

# For fun, but still useful
import smtplib
from email.mime.text import MIMEText
from clint.textui import puts, colored


def planetInString(planet):
    """
    Simbad does not like planet names, so remove 'b', 'c' etc. in the of the
    star name
    """
    letters = ['b', 'c', 'd', 'e', 'f', 'g']
    for letter in letters:
        if planet.endswith(' ' + letter):
            planet = planet[0:-2]
    return planet


def downloadExoplanet():
    link_csv = 'http://exoplanet.eu/catalog/csv/?f=%22controversial%22+IN+detection'

    puts(colored.clean('Updating from exoplanets.eu.'))
    puts(colored.yellow('Please wait...'))
    response = urllib2.urlopen(link_csv)
    table = response.read()
    with open("exo_cont.csv", "w") as f:
        f.write(table)

    exoplanet = csv.reader(open("exo_cont.csv"))
    exoplanet = sorted(exoplanet, key=operator.itemgetter(47), reverse=True)
    return exoplanet


def update(exoplanet, starsID, update_sweetcat):
    s = 0
    names = []
    for row in exoplanet[0:-1]:  # To avoid the header.
        tt = row[47]
        update_exoplanet = datetime(int(tt[0:4]), int(tt[5:7]), int(tt[8:10]))
        if update_exoplanet > update_sweetcat:
            updated_star = planetInString(row[0])
            if updated_star not in starsID:
                s += 1
                names.append(updated_star)
    return names


def sendingMail(names):
    # Sending the mail
    fp = open('mail.txt', 'rb')
    msg = MIMEText(fp.read())
    fp.close()

    with open('mailinfo.txt', 'rb') as f:
        sender = f.readline().split(': ')[1].strip('\n')
        receiver = f.readline().split(': ')[1].strip('\n')
        smtp = f.readline().split(': ')[1].strip('\n')

    msg['Subject'] = 'Update available to SWEET-Cat: ' + str(N) +\
        ' new controversial exoplanets'
    msg['From'] = sender
    msg['To'] = receiver

    try:
        s = smtplib.SMTP(smtp)
        s.sendmail(sender, [receiver], msg.as_string())
        s.quit()
    except smtplib.SMTPRecipientsRefused, e:
        puts(colored.red('Not able to send an email...'))


if __name__ == '__main__':
    SC = csv.reader(open('WEBSITE_online.rdb'), delimiter='\t')
    SC = sorted(SC, key=operator.itemgetter(24), reverse=True)
    starsID = []
    for i in xrange(len(SC)):
        starsID.append(SC[i][0])

    tmp = SC[0][-3]
    update_sweetcat = datetime(int(tmp[0:4]), int(tmp[5:7]), int(tmp[8:10]))

    exoplanet = downloadExoplanet()
    names = update(exoplanet, starsID, update_sweetcat)
    names = set(names)
    N = len(names)

    if N > 0:
        puts(colored.green(str(N) + " new exoplanet available!"))

        # Preparing list for SIMBAD
        simbad(names)

        sendingMail(names)
    else:
        puts(colored.clean('No new updates seems to be available.'))
        puts(colored.clean('SWEET-Cat should be up to date'))
