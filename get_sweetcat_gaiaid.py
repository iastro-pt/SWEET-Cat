#!/usr/bin/python
## My first python code

##imports:

import checkExoplanet as cEx
from astroquery.simbad import Simbad
import numpy as np
from astroquery.vizier import Vizier
from astropy import units as u

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

def create_sweetcat_gaiaid_rdb(filename = 'gaiadr2/gaiaid_sc.rdb'):
  customSimbad=custumized_Simbad()
  exodb = cEx.Update(controversial=False)
  print("Number of stars in SWEET-Cat: ", len(exodb.SC))

  data_gaia_manual = np.loadtxt("gaiadr2/gaia_lookup.rdb", skiprows=2,dtype="str",delimiter="\t").transpose()

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

def get_gaia_dr2_paralax(filename = 'gaiadr2/gaiaid_sc.rdb'):
  print("Result GAIA vizier")
  vq2 = Vizier(columns=['Source','Plx','e_Plx', 'FG','e_FG','Gmag','e_Gmag', 'BPmag','e_BPmag', 'RPmag','e_RPmag', 'o_Gmag'], row_limit=5000) 

  name, gaia_id = np.loadtxt(filename, unpack=True, usecols=(0,1), skiprows=2, delimiter="\t", dtype=str)
  radius_search = 10.0*u.arcsec

  strlines = ['name\tgaia_id\tPlx\te_Plx\tGmag\te_Gmag\tRPmag\te_RPmag\tBPmag\te_BPmag\tFG\te_FG\tG_flux_std_n\n']
  strlines.append('----\t-------\t---\t-----\t----\t------\t-----\t-------\t-----\t-------\t--\t----\t------------\n')
  for i,gaiadr2 in enumerate(gaia_id):
    print(i, len(gaia_id), name[i], gaiadr2)
    if gaiadr2 == "-1":
      print(name[i], "No gaia id dr2")
      strlines.append('%s\t%s\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\n' % (name[i], gaiadr2))
    else:
      result_gaia_vizier=vq2.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/345/gaia2"], radius=radius_search*15.)
      try:
        iline = np.where(result_gaia_vizier[0]['Source'] == int(gaiadr2))[0][0]
        print(result_gaia_vizier[0]['Source','Plx','e_Plx', 'Gmag', 'e_Gmag', 'FG', 'e_FG', 'RPmag', 'e_RPmag', 'BPmag', 'e_BPmag'][iline])
        std_g_flux_norm1_v = result_gaia_vizier[0]['e_FG'][iline]*np.sqrt(result_gaia_vizier[0]['o_Gmag'][iline])/result_gaia_vizier[0]['FG'][iline]
      except IndexError:
        result_gaia_vizier=vq2.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/345/gaia2"], radius=radius_search*25.)
        try:
          iline = np.where(result_gaia_vizier[0]['Source'] == int(gaiadr2))[0][0]
          print(result_gaia_vizier[0]['Source','Plx','e_Plx', 'Gmag', 'e_Gmag', 'FG', 'e_FG', 'RPmag', 'e_RPmag', 'BPmag', 'e_BPmag'][iline])
          std_g_flux_norm1_v = result_gaia_vizier[0]['e_FG'][iline]*np.sqrt(result_gaia_vizier[0]['o_Gmag'][iline])/result_gaia_vizier[0]['FG'][iline]
        except IndexError:
          result_gaia_vizier =  [Table( [[None], [None],[None], [None] ,[None] ,[None]   ,[None] , [None]    ] ,
                                   names=('Plx','e_Plx','Gmag','e_Gmag','RPmag','e_RPmag','BPmag','e_BPmag')) ]
          std_g_flux_norm1_v = -1
          iline=0
      print("Iline2:", iline)
      strlines.append('%s\t%s\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%.3e\t%.3e\t%10.7f\n' % (name[i], gaiadr2, 
        result_gaia_vizier[0]['Plx'][iline], result_gaia_vizier[0]['e_Plx'][iline],
        result_gaia_vizier[0]['Gmag'][iline], result_gaia_vizier[0]['e_Gmag'][iline],
        result_gaia_vizier[0]['RPmag'][iline], result_gaia_vizier[0]['e_RPmag'][iline],
        result_gaia_vizier[0]['BPmag'][iline], result_gaia_vizier[0]['e_BPmag'][iline],
        result_gaia_vizier[0]['FG'][iline], result_gaia_vizier[0]['e_FG'][iline],
        std_g_flux_norm1_v))
      print(strlines[-1])
  fileo = open('gaiadr2/tmp.rdb', "w")
  for line in strlines:
    fileo.write(line)
  fileo.close()


### Main program:
def main():
  print("Hello")
  #create_sweetcat_gaiaid_rdb()
  get_gaia_dr2_paralax()

if __name__ == "__main__":
    main()

