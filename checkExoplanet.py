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

SC = csv.reader(open('WEBSITE_online.rdb'), delimiter='\t')
SC = sorted(SC, key=operator.itemgetter(24), reverse=True)

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
for row in exoplanet:
    tt = row[38]
    if tt.startswith(' updated'):
        pass
    else:
        update_exoplanet = datetime(int(tt[0:4]), int(tt[5:7]), int(tt[8:10]))
        if update_exoplanet > update_sweetcat:
            s += 1
            names.append(row[0])

# Simbad does not like planet names, so remove 'b', 'c' etc. in the of the star
# name
letters = ['b', 'c', 'd', 'e', 'f', 'g']
for i in range(len(names)):
    for letter in letters:
        if names[i].endswith(' ' + letter):
            names[i] = names[i][0:-2]

if s > 0:
    print s, "new exoplanet available!"
    print "Sending mail to maintainer"

    # Preparing list for SIMBAD
    simbad(names)

    # Sending the mail
    fp = open('mail.txt', 'rb')
    msg = MIMEText(fp.read())
    fp.close()

    msg['Subject'] = 'Update available to SWEET-Cat: ' + str(s) +\
        ' new exoplanets'
    msg['From'] = 'daniel.andreasen@astro.up.pt'
    msg['To'] = 'daniel.andreasen@astro.up.pt'

    s = smtplib.SMTP('write smtp server here')
    s.sendmail('email from goes here',
               ['email to goes here'], msg.as_string())
    s.quit()
else:
    print "No new updates seems to be available."
    print "SWEET-Cat should be up to date"
