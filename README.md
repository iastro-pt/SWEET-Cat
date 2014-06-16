SWEET-Cat
=========

Repository for the SWEET-Cat(alog), a table for FGKM stars hosting exoplanets.

Here is the software for updating SWEET-Cat, and probably some of the latest
versions and the currently used version on the
[SWEET-Cat webpage](https://www.astro.up.pt/resources/sweet-cat/).
Note that the latest version is always available for download on the webpage.



Installation
============
It is now possible to get an update everytime there is a new planet on
[exoplanet.eu](http://www.exoplanet.eu/catalog).

To do this, place `SWEET-Cat` in your `PATH`. If you want an update everyday at
13 you can use crontab. Add the following line with `crontab -e`

    0 13 * * * /full/path/to/script/SWEET-Cat

Remember to change the receiver and sender in `checkExoplanet.py` line 45 and
46. You can also decide to change the `mail.txt` for whatever you like. This is
the text in the mail.

TODO
====

   - [ ] Make a function that look at the latest update on SWEET-Cat and compare
     with the latest update on exoplanet.eu. If not the same, SWEET-Cat should
     be updated.
   - Many of the SWEET-Cat things can be done automatically, e.g. write all the
     things in the table.
   - Create files automatically for ParralaxSpec.py, TorresMass.py, and
     Simbad.py.
   - Create a NEWNEW1. This one is just a file with the names for creating the
     Simbad script (see e.g. in folder).
   - The toughest challenge must be to find values in articles. I guess we must
     do this part manually, but we might be able to make a script that notifies
     us, if with finds article from Santos or Sousa, so we can set the homogenity
     flag to 1.
