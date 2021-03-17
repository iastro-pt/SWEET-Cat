#!/usr/bin/python
## My first python code

##imports:
import sys
sys.path.insert(0, '../') 
import checkExoplanet as cEx
from astroquery.simbad import Simbad
import numpy as np
from astroquery.vizier import Vizier
from astropy import units as u
import os
from astropy.table import Table
import logg
import math

## My functions:
def custumized_Simbad():
  customSimbad = Simbad()
  customSimbad.add_votable_fields('flux(G)') 
  customSimbad.add_votable_fields('flux(B)') 
  customSimbad.add_votable_fields('flux(V)') 
  customSimbad.add_votable_fields('flux(R)') 
  customSimbad.add_votable_fields('flux(I)') 
  customSimbad.add_votable_fields('flux(J)') 
  customSimbad.add_votable_fields('flux(H)') 
  customSimbad.add_votable_fields('flux(K)')
  customSimbad.add_votable_fields('flux_error(G)') 
  customSimbad.add_votable_fields('flux_error(B)') 
  customSimbad.add_votable_fields('flux_error(V)') 
  customSimbad.add_votable_fields('flux_error(R)') 
  customSimbad.add_votable_fields('flux_error(I)') 
  customSimbad.add_votable_fields('flux_error(J)') 
  customSimbad.add_votable_fields('flux_error(H)') 
  customSimbad.add_votable_fields('flux_error(K)')
  customSimbad.add_votable_fields('plx')
  customSimbad.add_votable_fields('plx_error')
  return customSimbad

def get_gaia_dr2_id(results_ids):
  for name in results_ids[::-1]:
    #print(name[0])
    if "Gaia DR2 " in name[0]:
      return name[0].split(" ")[-1]
  return -1

def create_sweetcat_gaiaid_rdb(filename = 'gaiaid_sc.rdb'):
  customSimbad=custumized_Simbad()
  exodb = cEx.Update(controversial=False)
  print("Number of stars in SWEET-Cat: ", len(exodb.SC))

  data_gaia_manual = np.loadtxt("gaia_lookup.rdb", skiprows=2,dtype="str",delimiter="\t").transpose()

  strlines = ['name\tgaia_id\n','----\t-------\n']

  for i in range(len(exodb.SC)):
    name = exodb.SC.iloc[i]['name']
    name_o = name
    print(i, name)
    if name[-2:] == " A":
        name =  name[:-2]
    if "(AB)" in name:
        name = name.replace("(AB)", "")
    if "Qatar" in name:
        name = name.replace("-","")
    try:
        result_ids = customSimbad.query_objectids(name)
    except:
        result_ids = customSimbad.query_objectids(name)
    if result_ids is None:
        gaiadr2 = -1
    else:
        gaiadr2 = get_gaia_dr2_id(result_ids)
    if gaiadr2 == -1:
        idgm = np.where(name_o == data_gaia_manual[0])[0][0]
        gaiadr2 = data_gaia_manual[2][idgm]
    print(gaiadr2)
    strlines.append('%s\t%s\n' % (name_o, gaiadr2))

  fileo = open(filename, "w")
  for line in strlines:
    fileo.write(line)
  fileo.close()

def get_gaia_dr3_paralax(filename = 'gaiaid_sc.rdb', fileout= 'sweet_gaia_dr3.rdb', append= True):
  if os.path.isfile(fileout):
    filetmp = open(fileout,"r")
    strlines = filetmp.readlines()
    filetmp.close()
    name_tmp, gaiaid_tmp = np.loadtxt(fileout, unpack=True,usecols=(0,1), skiprows=2, delimiter="\t", dtype=str)
  else:
    print("Fresh start")
    gaiaid_tmp = []
    strlines =     ['name\tgaia_id\tPlx\te_Plx\tGmag\te_Gmag\tRPmag\te_RPmag\tBPmag\te_BPmag\tFG\te_FG\tG_flux_std_n\n']
    strlines.append('----\t-------\t---\t-----\t----\t------\t-----\t-------\t-----\t-------\t--\t----\t------------\n')

  print("Result GAIA vizier")
  vq2 = Vizier(columns=['Source','Plx','e_Plx', 'FG','e_FG','Gmag','e_Gmag', 'BPmag','e_BPmag', 'RPmag','e_RPmag', 'o_Gmag'], row_limit=5000) 

  name, gaia_id = np.loadtxt(filename, unpack=True, usecols=(0,1), skiprows=2, delimiter="\t", dtype=str)
  radius_search = 10.0*u.arcsec

  #ist = 3100
  #for i,gaiadr2 in enumerate(gaia_id[ist:]):
  #  i += ist
  for i,gaiadr2 in enumerate(gaia_id):
    print(i, len(gaia_id), name[i], gaiadr2)
    if gaiadr2 in gaiaid_tmp and name[i] in name_tmp:
      print("This one already in place")
      continue
    if gaiadr2 == "-1":
      print(name[i], "No gaia id dr2")
      strlines.append('%s\t%s\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\n' % (name[i], gaiadr2))
    else:
      result_gaia_vizier_dr3=vq2.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/350/gaiaedr3"], radius=radius_search*15.)
      try:
        iline3 = np.where(result_gaia_vizier_dr3[0]['Source'] == int(gaiadr2))[0][0]
        print(result_gaia_vizier_dr3[0]['Source','Plx','e_Plx', 'Gmag', 'e_Gmag', 'FG', 'e_FG', 'RPmag', 'e_RPmag', 'BPmag', 'e_BPmag'][iline3])
        std_g_flux_norm1_v = result_gaia_vizier_dr3[0]['e_FG'][iline3]*np.sqrt(result_gaia_vizier_dr3[0]['o_Gmag'][iline3])/result_gaia_vizier_dr3[0]['FG'][iline3]
      except IndexError:
        result_gaia_vizier_dr3=vq2.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/350/gaiaedr3"], radius=radius_search*25.)
        try:
          iline3 = np.where(result_gaia_vizier_dr3[0]['Source'] == int(gaiadr2))[0][0]
          print(result_gaia_vizier_dr3[0]['Source','Plx','e_Plx', 'Gmag', 'e_Gmag', 'FG', 'e_FG', 'RPmag', 'e_RPmag', 'BPmag', 'e_BPmag'][iline3])
          std_g_flux_norm1_v = result_gaia_vizier_dr3[0]['e_FG'][iline3]*np.sqrt(result_gaia_vizier_dr3[0]['o_Gmag'][iline3])/result_gaia_vizier_dr3[0]['FG'][iline3]

        except IndexError:
          result_gaia_vizier_dr3 =  [Table( [[-1], [-1],[-1], [-1] ,[-1], [-1] ,[-1] ,[-1]   ,[-1] , [-1]    ] ,
                                   names=('Plx','e_Plx','Gmag','e_Gmag','FG','e_FG','RPmag','e_RPmag','BPmag','e_BPmag')) ]
          std_g_flux_norm1_v = -1
          iline3=0
      print("Iline2:", iline3)
      strlines.append('%s\t%s\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%.3e\t%.3e\t%10.7f\n' % (name[i], gaiadr2, 
        result_gaia_vizier_dr3[0]['Plx'][iline3], result_gaia_vizier_dr3[0]['e_Plx'][iline3],
        result_gaia_vizier_dr3[0]['Gmag'][iline3], result_gaia_vizier_dr3[0]['e_Gmag'][iline3],
        result_gaia_vizier_dr3[0]['RPmag'][iline3], result_gaia_vizier_dr3[0]['e_RPmag'][iline3],
        result_gaia_vizier_dr3[0]['BPmag'][iline3], result_gaia_vizier_dr3[0]['e_BPmag'][iline3],
        result_gaia_vizier_dr3[0]['FG'][iline3], result_gaia_vizier_dr3[0]['e_FG'][iline3],
        std_g_flux_norm1_v))
      print(strlines[-1])
  if append:     
    fileo = open(fileout, "a")
  else:
    fileo = open(fileout, "w")
  for line in strlines:
    fileo.write(line)
  fileo.close()

def get_logg_gaia(filename = 'sweet_gaia_dr3.rdb', fileout='sweet_logg_dr3.rdb', append= True):
  if os.path.isfile(fileout):
    filetmp = open(fileout,"r")
    strlines = filetmp.readlines()
    filetmp.close()
    name_tmp, gaiaid_tmp = np.loadtxt(fileout, unpack=True,usecols=(0,1), skiprows=2, delimiter="\t", dtype=str)
  else:
    print("Fresh start")
    gaiaid_tmp = []
    strlines =     ['name\tgaia_id\tlogg_gaia\te_logg_gaia\tlogg_sw\te_logg_sw\tteff_sw\te_teff_sw\tsource\n']
    strlines.append('----\t-------\t---------\t-----------\t-------\t---------\t-------\t---------\t------\n')

  name, gaia_id = np.loadtxt(filename, unpack=True, usecols=(0,1), skiprows=2, delimiter="\t", dtype=str)
  plx,eplx, gmag, egmag = np.loadtxt(filename, unpack=True, usecols=(2,3,4,5), skiprows=2, delimiter="\t")
  exodb = cEx.Update(controversial=False)

  for i,gaiadr2 in enumerate(gaia_id):
    mass = exodb.SC.iloc[i]['M']
    ermass = exodb.SC.iloc[i]['Merr']
    teff = exodb.SC.iloc[i]['Teff']
    erteff = exodb.SC.iloc[i]['Tefferr']
    loggsw = exodb.SC.iloc[i]['logg']
    erloggsw = exodb.SC.iloc[i]['logger']
    source = exodb.SC.iloc[i]['source']
    if math.isnan(erteff):
      print("error on Teff None")
      erteff = 250

    print(i, len(gaia_id), name[i], gaiadr2, gmag[i], plx[i], exodb.SC.iloc[i]['name'], teff, erteff, mass, ermass)
    if gaiadr2 in gaiaid_tmp and name[i] in name_tmp:
      print("This one already in place")
      continue
    if gaiadr2 == "-1" or gmag[i] == -1 or math.isnan(mass) or math.isnan(teff):
      if math.isnan(mass) or math.isnan(teff):
        print("no mass or teff in sweet-cat")
        teff=-1; erteff=-1
      else:
        print(name[i], "No gaia id dr2")
      strlines.append('%s\t%s\t-1\t-1\t%5.2f\t%5.2f\t%6d\t%6d\t%2d\n' % (name[i], gaiadr2, loggsw, erloggsw, teff, erteff, source))
    else:
      gaia_gmag = gmag[i]
      ergaia_gmag = egmag[i]
      gaia_paralax = plx[i]
      ergaia_paralax = eplx[i]
      print(mass, ermass, teff, erteff, gaia_gmag, ergaia_gmag, gaia_paralax, ergaia_paralax)
      logg_gaiadr3, erlogg_gaiadr3 = logg.logg_gaia_error(mass, ermass, teff, erteff, gaia_gmag, ergaia_gmag, gaia_paralax, ergaia_paralax)
      print(logg_gaiadr3, erlogg_gaiadr3, loggsw, erloggsw)
      strlines.append('%s\t%s\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%6d\t%6d\t%2d\n' % (name[i], gaiadr2, logg_gaiadr3, erlogg_gaiadr3, loggsw, erloggsw, teff, erteff, source))
  fileo = open(fileout, "w")
  for line in strlines:
    fileo.write(line)
  fileo.close()






### Main program:
def main():
  print("Hello")
  #create_sweetcat_gaiaid_rdb()
  #get_gaia_dr3_paralax()
  get_logg_gaia()


if __name__ == "__main__":
    main()

