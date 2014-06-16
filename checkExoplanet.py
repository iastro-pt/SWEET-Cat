#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
import csv
import operator
import urllib2
from datetime import datetime

# For fun, but still useful
import smtplib
from email.mime.text import MIMEText

# Link to download the latest table from exoplanets.eu
link_csv = 'http://www.exoplanet.eu/catalog/csv/?f=%22radial%22+IN+detection'
link_csv += '+OR+%22astrometry%22+IN+detection+OR+%22transit%22+IN+detection'

print "Updating from exoplanets.eu.\nPlease wait..."
response = urllib2.urlopen(link_csv)
table = response.read()
with open("exo.csv", "w") as f:
    f.write(table)
exoplanet = csv.reader(open("exo.csv"))
exoplanet = sorted(exoplanet, key=operator.itemgetter(38), reverse=True)


tmp = exoplanet[0][38]
update_exoplanet = datetime(int(tmp[0:4]), int(tmp[5:7]), int(tmp[8:10]))

SC = csv.reader(open('WEBSITE_online.rdb'), delimiter='\t')
SC = sorted(SC, key=operator.itemgetter(24), reverse=True)

tmp = SC[0][-3]
update_sweetcat = datetime(int(tmp[0:4]), int(tmp[5:7]), int(tmp[8:10]))

if update_exoplanet > update_sweetcat:
    print "New exoplanet available!"
    print "Sending mail to maintainer"

    fp = open('mail.txt', 'rb')
    msg = MIMEText(fp.read())
    fp.close()

    msg['Subject'] = 'Update available to SWEET-Cat'
    msg['From'] = 'daniel.andreasen@astro.up.pt'
    msg['To'] = 'daniel.andreasen@astro.up.pt'

    s = smtplib.SMTP('mail.astro.up.pt')
    s.sendmail('daniel.andreasen@astro.up.pt',
               ['daniel.andreasen@astro.up.pt'], msg.as_string())
    s.quit()
else:
    print "No new updates seems to be available."
    print "SWEET-Cat should be up to date"
