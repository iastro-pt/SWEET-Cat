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


SC = csv.reader(open('WEBSITE_online.rdb'), delimiter='\t')
SC = sorted(SC, key=operator.itemgetter(24), reverse=True)
starsID = []
for i in range(len(SC)):
    starsID.append(SC[i][0])


# Link to download the latest table from exoplanets.eu
tmp = SC[0][-3]
update_sweetcat = datetime(int(tmp[0:4]), int(tmp[5:7]), int(tmp[8:10]))
link_csv = 'http://www.exoplanet.eu/catalog/csv/?f=%22radial%22+IN+detection'
link_csv += '+OR+%22astrometry%22+IN+detection+OR+%22transit%22+IN+detection'

print "Updating from exoplanets.eu.\nPlease wait..."
response = urllib2.urlopen(link_csv)
table = response.read()
with open("exo.csv", "w") as f:
    f.write(table)
exoplanet = csv.reader(open("exo.csv"))
exoplanet = sorted(exoplanet, key=operator.itemgetter(38), reverse=True)


s = 0
names = []
for row in exoplanet[0:-1]:  # To avoid the header.
    tt = row[38]
    update_exoplanet = datetime(int(tt[0:4]), int(tt[5:7]), int(tt[8:10]))
    if update_exoplanet > update_sweetcat:
        updated_star = planetInString(row[0])
        if updated_star not in starsID:
            s += 1
            names.append(updated_star)


if s > 0:
    print s, "new exoplanet available!"
    print "Sending mail to maintainer"

    # Preparing list for SIMBAD
    simbad(set(names))

    # Sending the mail
    fp = open('mail.txt', 'rb')
    msg = MIMEText(fp.read())
    fp.close()

    with open('mailinfo.txt', 'rb') as f:
        sender = f.readline().split(': ').strip('\n')
        receiver = f.readline().split(': ').strip('\n')
        smtp = f.readline().split(': ').strip('\n')

    msg['Subject'] = 'Update available to SWEET-Cat: ' + str(s) +\
        ' new exoplanets'
    msg['From'] = sender
    msg['To'] = receiver

    s = smtplib.SMTP(smtp)
    s.sendmail(sender, [receiver], msg.as_string())
    s.quit()
else:
    print "No new updates seems to be available."
    print "SWEET-Cat should be up to date"
