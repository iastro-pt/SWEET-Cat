SWEET-Cat
=========

Repository for the SWEET-Cat(alog), a table for FGKM stars hosting exoplanets.

Here is the software for updating SWEET-Cat, and probably some of the latest
versions and the currently used version on the
[SWEET-Cat webpage](https://www.astro.up.pt/resources/sweet-cat/).
Note that the latest version is always available for download on the web page.



Installation
============
It is now possible to get an update every time there is a new planet on
[exoplanet.eu](http://www.exoplanet.eu/catalog).

To do this, place `SWEET-Cat` in your `PATH`. If you want an update everyday at
13 you can use crontab. Add the following line with `crontab -e`

    0 13 * * * /full/path/to/script/SWEET-Cat

Remember to change the receiver and sender in `checkExoplanet.py` line 45 and
46. You can also decide to change the `mail.txt` for whatever you like. This is
the text in the mail.

Setting up the mail
===================
A file names `mailinfo.txt` needs to be created. It should looks like the
following

    sender: from@host.com
    receiver: to@hotst.com
    smtp: mail.com

Change the three options with what works for you. Questions can be directed to
the author.


TODO
====

   - [x] Make a function that look at the latest update on SWEET-Cat and compare
     with the latest update on exoplanet.eu. If not the same, SWEET-Cat should
     be updated.
   - [x] Find how many new exoplanets is found and send that in email as well.
   - [x] Make a good interface between `checkExoplanet.py` (which seems to be the
     main function) and the other scripts. `os.system("python Simbad.py
     `NEWNEW1`")` is not nice...
   - [x] Create a `NEWNEW1`. This one is just a file with the names for creating the
     Simbad script (see e.g. in folder).
   - [x] It is also important to check for false-positives, so we can remove
     them as well.
   - [x] When updating, make sure that the star is not already in SWEET-Cat. If
     that is the case, just look at what have been updated.
   - [x] If the star is already in SWEET-Cat, find a way to deal with it so I
     don't get notifications all the time...
   - [x] Would be nice if the two scripts for checking for updates have a
     common file with information necessary for sending the mail, like the
     sender, reciever, and smpt information.
   - [x] Remove duplicate planet host, if e.g. two exoplanets have been
     detected around the same host, the host should only appear once in the
     `NEWNEW1` file.
   - [ ] Many of the SWEET-Cat things can be done automatically, e.g. write all the
     things in the table.
   - [ ] Create files automatically for `ParralaxSpec.py`, `TorresMass.py`, and
     Simbad.py.
   - [ ] The toughest challenge must be to find values in articles. I guess we must
     do this part manually, but we might be able to make a script that notifies
     us, if with finds article from Santos or Sousa, so we can set the homogenity
     flag to 1.

Stuff added by Sérgio: (not necessarily to appear online for now)

   - [ ] Include a column in the table with a flag (yes or no) to inform that a
     given planet-host has transit events detected
   - [ ] The interface that we use to get the mass: stev.oapd.inaf.it/cgi-bin/param_1.3
     also provides radius and ages. It may be useful to include this in the table
