#!/usr/bin/env python
# -*- coding: utf-8 -*-

# My imports
import numpy as np
from astropy import coordinates as coord
from astropy import units as u
import pandas as pd
from clint.textui import puts, colored
import time
from astroquery.simbad import Simbad
import warnings
warnings.filterwarnings('ignore')
from astroquery.irsa_dust import IrsaDust
from astroquery.vizier import Vizier
from astropy.table import Table
import math

import sys
#sys.path.insert(0, "/home/sousasag/Programas/GIT_Projects/MR_spec/")
sys.path.insert(0, "/home/sousasag/Programas/GIT_projects/MR_spec/")
import logg_paralax_mass as loggdr3
import MR_spec as MR




columns_SW = ['Name', 'hd', 'RA', 'DEC', 'Vmag', 'eVmag', 'PlxFlag', 'Teff', 'eTeff',
       'Logg', 'eLogg', 'Vt', 'eVt', '[Fe/H]', 'e[Fe/H]', 'Reference', 'Link',
       'SWFlag', 'Update', 'Comment', 'Database', 'gaia_dr2', 'gaia_dr3',
       'Plx', 'ePlx', 'Gmag', 'eGmag', 'RPmag', 'eRPmag', 'BPmag', 'eBPmag',
       'FG', 'eFG', 'G_flux_std_n', 'Logg_gaia', 'eLogg_gaia', 'Mass_t',
       'eMass_t', 'Radius_t', 'eRadius_t', 'spec_base', 'Distance', 'RA_EU',
       'DEC_EU', 'RA_NASA', 'DEC_NASA', 'Distance_b', 'eDistance_b']


dtype_SW = {'gaia_dr2':'int64','gaia_dr3':'int64'}

def remove_planet(name):
    """Remove the trailing b, c, d, etc in the stellar name"""
    planets = 'abcdefghijB'
    for planet in planets:
        if name.endswith(' %s' % planet):
            return name[:-2].strip()
    # some exoplanets have .01 or .02 in the name 
    if name.endswith('.01') or name.endswith('.02') or name.endswith('.2'):
        return name[:-3].strip()    
    if  name.endswith(' (bc)'):
        return name[:-4].strip()
    return name

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

def get_gaiadr2(name):
    customSimbad=custumized_Simbad()
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
    return gaiadr2

def get_gaia_dr3_id(results_ids):
  for name in results_ids[::-1]:
    #print(name[0])
    if "Gaia DR3 " in name[0]:
      return name[0].split(" ")[-1]
  return -1

def get_gaiadr3(name):
    customSimbad=custumized_Simbad()
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
        gaiadr3 = -1
    else:
        gaiadr3 = get_gaia_dr3_id(result_ids)
    return gaiadr3


def get_gaiadr3_data(gaiadr2, gaia3, name_search=None):
    vq2 = Vizier(columns=['Source','Plx','e_Plx', 'FG','e_FG','Gmag','e_Gmag', 'BPmag','e_BPmag', 'RPmag','e_RPmag', 'o_Gmag'], row_limit=5000) 
    radius_search = 10.0*u.arcsec
    if name_search is None:
        result_gaia_vizier_dr3=vq2.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/350/gaiaedr3"], radius=radius_search*15.)
    else:
        result_gaia_vizier_dr3=vq2.query_object(name_search, catalog=["I/350/gaiaedr3"], radius=radius_search*15.)
    try:
        iline3 = np.where(result_gaia_vizier_dr3[0]['Source'] == int(gaia3))[0][0]
        print(result_gaia_vizier_dr3[0]['Source','Plx','e_Plx', 'Gmag', 'e_Gmag', 'FG', 'e_FG', 'RPmag', 'e_RPmag', 'BPmag', 'e_BPmag'][iline3])
        std_g_flux_norm1_v = result_gaia_vizier_dr3[0]['e_FG'][iline3]*np.sqrt(result_gaia_vizier_dr3[0]['o_Gmag'][iline3])/result_gaia_vizier_dr3[0]['FG'][iline3]
    except IndexError:
        result_gaia_vizier_dr3=vq2.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/350/gaiaedr3"], radius=radius_search*25.)
        try:
            iline3 = np.where(result_gaia_vizier_dr3[0]['Source'] == int(gaia3))[0][0]
            print(result_gaia_vizier_dr3[0]['Source','Plx','e_Plx', 'Gmag', 'e_Gmag', 'FG', 'e_FG', 'RPmag', 'e_RPmag', 'BPmag', 'e_BPmag'][iline3])
            std_g_flux_norm1_v = result_gaia_vizier_dr3[0]['e_FG'][iline3]*np.sqrt(result_gaia_vizier_dr3[0]['o_Gmag'][iline3])/result_gaia_vizier_dr3[0]['FG'][iline3]

        except IndexError:
            result_gaia_vizier_dr3 =  Table( [[-1], [-1],[-1], [-1] ,[-1], [-1] ,[-1] ,[-1]   ,[-1] , [-1]    ] ,
                                   names=('Plx','e_Plx','Gmag','e_Gmag','FG','e_FG','RPmag','e_RPmag','BPmag','e_BPmag'))
            return result_gaia_vizier_dr3, -1
    print("Iline2:", iline3)
    print(result_gaia_vizier_dr3)
    return result_gaia_vizier_dr3[0][iline3], std_g_flux_norm1_v

def get_gaiadr3_dist(gaiadr2, gaia3):
    #vqdr3 = Vizier(columns=['Source','Plx','e_Plx', 'FG','e_FG','Gmag','e_Gmag', 'BPmag','e_BPmag', 'RPmag','e_RPmag', 'o_Gmag'], row_limit=5000) 
    vqdr3_dist = Vizier(row_limit=5000)
    radius_search = 10.0*u.arcsec
    try:
      #result_gaia_vizier_dr3=vqdr3.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/350/gaiaedr3"], radius=radius_search)
      result_dist = vqdr3_dist.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/352/gedr3dis"], radius=radius_search)
      #iline = np.where(result_gaia_vizier_dr3[0]['Source'] == int(gaia3))[0][0]
      #print(result_gaia_vizier_dr3[0]['Source','Plx','e_Plx', 'Gmag', 'e_Gmag', 'FG', 'e_FG', 'RPmag', 'e_RPmag', 'BPmag', 'e_BPmag'][iline])
#    std_g_flux_norm1_v = result_gaia_vizier_dr3[0]['e_FG'][iline]*np.sqrt(result_gaia_vizier_dr3[0]['o_Gmag'][iline])/result_gaia_vizier_dr3[0]['FG'][iline]
      iline_d = np.where(result_dist[0]['Source'] == int(gaia3))[0][0]
      #print(result_dist[0]['Source','rgeo', 'rpgeo'][iline_d])
    except:
      try:
        #result_gaia_vizier_dr3=vqdr3.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/350/gaiaedr3"], radius=100*radius_search)
        result_dist = vqdr3_dist.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/352/gedr3dis"], radius=20*radius_search)
        #iline = np.where(result_gaia_vizier_dr3[0]['Source'] == int(gaia3))[0][0]
        iline_d = np.where(result_dist[0]['Source'] == int(gaia3))[0][0]
      except:
        print("problem")
        return -1,-1,-1,-1
    #print("Ilines:", iline, iline_d)
    rgeo = result_dist[0][iline_d]["rgeo"]
    b_rgeo = result_dist[0][iline_d]["b_rgeo"]
    B_rgeo = result_dist[0][iline_d]["B_rgeo"]
    ergeo = np.max([rgeo-b_rgeo, B_rgeo-rgeo])
    rpgeo = result_dist[0][iline_d]["rgeo"]
    b_rpgeo = result_dist[0][iline_d]["b_rpgeo"]
    B_rpgeo = result_dist[0][iline_d]["B_rpgeo"]
    erpgeo = np.max([rpgeo-b_rpgeo, B_rpgeo-rpgeo])
    return rgeo, ergeo, rpgeo, erpgeo



def get_logg_dr3_stuff(Teff, Tefferr, logg, loggerr, FeH, FeHerr, Gmag, e_Gmag, Plx, e_Plx):
    try:
        Teff = float(Teff)
    except:
        print("No Temperature")
        return -1, -1, -1, -1, -1, -1
    Tefferr = 200 if Tefferr == "NULL" else float(Tefferr)
    if math.isnan(Tefferr):
        Tefferr = 200
    FeH = 0 if FeH == "NULL" else float(FeH) 
    if math.isnan(FeH):
        FeH = 0
    FeHerr = 0.2 if FeHerr == "NULL" else float(FeHerr)
    if math.isnan(FeHerr):
        FeHerr=0.2
    logg = 4.4 if logg == "NULL" else float(logg)    
    if math.isnan(logg):
        logg = 4.4
    loggerr = 0.5 if loggerr == "NULL" else float(loggerr)
    if math.isnan(loggerr):
        loggerr=0.5
    try:
        logg_gaiadr3, erlogg_gaiadr3, masst, emasst = loggdr3.logg_mass_iteractive_error(Teff, Tefferr, logg, loggerr, FeH, FeHerr, Gmag, e_Gmag, Plx, e_Plx)
        print(Teff, Tefferr, logg, loggerr, FeH, FeHerr, Gmag, e_Gmag, Plx, e_Plx)
        radiust, eradiust = MR.radius_torres2010_error(Teff, Tefferr, logg_gaiadr3, erlogg_gaiadr3, FeH, FeHerr)
    except:
        return -1, -1, -1, -1, -1, -1
    return logg_gaiadr3, erlogg_gaiadr3, masst, emasst, radiust, eradiust






##Parallax Spec:
def bolcor(teff):
    """
    Calculate the bolometric correction, given the temperature
    """
    lteff = np.log10(teff)
    bcflow=np.nan

    if lteff < 3.7:
        bcflow = -0.190537291496456e+05 + 0.155144866764412e+05 * lteff -\
            0.421278819301717e+04 * (lteff * lteff) +\
            0.381476328422343e+03 * (lteff*lteff*lteff)

    if lteff >= 3.7 and lteff < 3.9:
        bcflow = -0.370510203809015e+05 + 0.385672629965804e+05 * lteff -\
            0.150651486316025e+05 * (lteff * lteff) +\
            0.261724637119416e+04 * (lteff*lteff*lteff) -\
            0.170623810323864e+03 * (lteff * lteff * lteff * lteff)
    else:
        bcflow = -0.118115450538963e+06 + 0.137145973583929e+06 * lteff -\
            0.636233812100225e+05 * (lteff * lteff) +\
            0.147412923562646e+05 * (lteff * lteff * lteff) -\
            0.170587278406872e+04 * (lteff * lteff * lteff * lteff) +\
            0.788731721804990e+02 * (lteff * lteff * lteff * lteff * lteff)
    
    return bcflow


def parallax(teff,eteff, logg,elogg,vmag,evmag,   mass,emass,  Av,eAv):
    """
    Calculate the parallax, given the mass Santos 2004
    """
    if eteff=='NULL' or elogg=='NULL' or evmag=='NULL' or  emass=='NULL' or  eAv=='NULL':
        bcflow = bolcor(teff)
        return 10.**((logg - 4.44 - np.log(mass)/np.log(10.) - 4.*np.log(teff) / np.log(10.) + \
            4.*np.log(5777.)/np.log(10.) - 0.4*(vmag + bcflow-Av) - 0.11) * 0.5) * 1000, np.nan

    ntrials = 10000
    randomteff = float(teff) + float(eteff)*np.random.randn(ntrials)
    randomlogg = float(logg) + float(elogg)*np.random.randn(ntrials)
    randommass = abs(mass + emass*np.random.randn(ntrials))
    randomvmag = float(vmag) + float(evmag)*np.random.randn(ntrials)
    randomAv = Av + eAv*np.random.randn(ntrials)
    par=[]
    for i in range(ntrials):
        bcflow = bolcor(randomteff[i])
        par.append(10.**((randomlogg[i] - 4.44 - np.log(randommass[i])/np.log(10.) \
        - 4.*np.log(randomteff[i]) / np.log(10.) + 4.*np.log(5777.)/np.log(10.) - 0.4*\
        (randomvmag[i] + bcflow-randomAv[i]) - 0.11) * 0.5) * 1000)
    
    sig=np.sqrt(np.sum((par - np.mean(par))**2) / (ntrials - 1))
    return np.mean(par),sig    





def GAIAplx(ra, de):
    v = Vizier(columns=["*", "+_r"], catalog='I/345/gaia2')
    pos=coord.SkyCoord(ra=ra, dec=de,unit=(u.hourangle,u.deg),frame='icrs',obstime='J2000')
    result=v.query_region(pos, radius="10s", catalog='I/345/gaia2')
    # Moving the positions to 2000
    try:
        nlines=len(result[0]['RA_ICRS'])
        deltat=-15.5
        sep=[]
        for ig,name in enumerate(result[0]['Source']):
            raold=result[0]['RA_ICRS'].data[ig]+(result[0]['pmRA'].data[ig] *deltat)/3600000.
            deold=result[0]['DE_ICRS'].data[ig]+(result[0]['pmDE'].data[ig] *deltat)/3600000.
            posold = coord.ICRS(ra=raold * u.deg, dec=deold * u.deg)
            sep.append(pos.separation(posold).arcsecond)
        indG=np.argmin(sep)
        if sep[indG]<1.5 and result[0]['Plx'].data[indG]>0:
            return str(round(result[0]['Plx'].data[indG],2)), str(round(result[0]['e_Plx'].data[indG],2))
    except:
        return 'NULL','NULL'

    return 'NULL','NULL'




def update_row(old_row, gaia_dr2, gaia_dr3, name_search = None):
    new_row = old_row.copy()
    new_row["gaia_dr2"] = gaia_dr2
    new_row["gaia_dr3"] = gaia_dr3
    gaia_data, std_g_flux = get_gaiadr3_data(gaia_dr2, gaia_dr3, name_search=name_search)
    rgeo, ergeo, rpgeo, erpgeo = get_gaiadr3_dist(gaia_dr2, gaia_dr3)

    teff     = old_row.Teff
    erteff   = old_row.eTeff
    feh      = old_row["[Fe/H]"]
    erfeh    = old_row["e[Fe/H]"]
    loggsw   = old_row["Logg"]
    erloggsw = old_row["eLogg"]
    gmag     = gaia_data['Gmag']
    egmag    = gaia_data['e_Gmag']
    gaia_paralax  = gaia_data['Plx']
    egaia_paralax = gaia_data['e_Plx']

    if math.isnan(erteff):
        print("error on Teff None")
        erteff = 250   
    if math.isnan(feh):
        feh = 0
    if math.isnan(erfeh):
        erfeh=0.2
    if math.isnan(loggsw):
        loggsw = 4.4
    if math.isnan(erloggsw):
        erloggsw=0.5
    meanlogg, stdlogg, meanmass, stdmass, r, er = loggdr3.get_logg_mass_radius_gaia_torres(teff, erteff, loggsw, erloggsw, feh, erfeh, gmag, egmag, gaia_paralax, egaia_paralax, distance_gaia=rgeo, e_distance_gaia=ergeo)

    new_row["PlxFlag"]      = "GAIAeDR3"
    new_row["gaia_dr2"]     = gaia_dr2
    new_row["gaia_dr3"]     = gaia_dr3

    new_row["Plx"]          = gaia_data["Plx"]
    new_row["e_Plx"]        = gaia_data["e_Plx"] 
    new_row["Gmag"]         = gaia_data["Gmag"]
    new_row["e_Gmag"]       = gaia_data["e_Gmag"]
    new_row["RPmag"]        = gaia_data["RPmag"]
    new_row["e_RPmag"]      = gaia_data["e_RPmag"] 
    new_row["BPmag"]        = gaia_data["BPmag"]
    new_row["e_BPmag"]      = gaia_data["e_BPmag"]
    new_row["FG"]           = gaia_data["FG"]
    new_row["e_FG"]         = gaia_data["e_FG"]
    new_row["G_flux_std_n"] = std_g_flux

    new_row["logg_gaia"]    = meanlogg
    new_row["e_logg_gaia"]  = stdlogg
    new_row["mass_t"]       = meanmass
    new_row["e_mass_t"]     = stdmass
    new_row["radius_t"]     = r
    new_row["e_radius_t"]   = er

    new_row["Distance"]     = 1000./new_row["Plx"].astype(float)

    new_row["Distance_b"]   = rgeo
    new_row["eDistance_b"]  = ergeo

    return new_row



def variable_assignment(digits):
    try:
        if digits > 0:
            x = '%.2f' % round(eval(input('> ')), digits)
        else:
            x = '%d' % round(eval(input('> ')), digits)
    except SyntaxError as e:
        x = 'NULL'
    return x

def read_float():
    val = input('> ')
    try:
        val = float(val)
    except:
        val = -1
    return val


if __name__ == '__main__':

    nasa = False
    nasa_str = input('\nExoplanet from NASA database? [True/False]: ')
    print(type(nasa_str), nasa_str)
    if nasa_str == 'True':
        nasa = True
        print('Adding exoplanets from the NASA database')
    elif nasa_str == 'False':
        nasa = False
        print('Adding exoplanets from the EU database')
    else:
        print('Answer different from True or False\nBye...')
        exit

    # stars = np.loadtxt('names.txt', dtype='S', delimiter='\t',usecols=(0,), )	"does not work when the file has only one line"
#    with open('names.txt') as f:
#        stars = f.readlines()
#    f.close()
    print(nasa)
    filenamesin = "names.txt"
    if nasa:
        filenamesin = "namesnasa.txt"
        #filenamesin = "namesnasa_c.txt"
    filenames = open(filenamesin, "r")
    linesnames = filenames.readlines()
    filenames.close()

    stars = np.loadtxt(filenamesin, delimiter='\t',usecols=(0,), dtype=str)
    ra_stars, dec_stars = np.loadtxt(filenamesin, delimiter='\t',usecols=(1,2), unpack=True)
    #print(stars, ra_stars, dec_stars)
    var = 'Y'

    # Read the data from exoplanet.eu
    fields = ['name', 'star_name', 'ra', 'dec', 'mag_v',
              'star_metallicity', 'star_metallicity_error_min',
              'star_metallicity_error_max',
              'star_teff', 'star_teff_error_min', 'star_teff_error_max']

    # Read data from NASA exoplanet archive
    fields_nasa = ['pl_hostname', 'hd_name', 'ra', 'dec', 'ra_str', 'dec_str',
                   'st_vj', 'st_vjerr',
                   'st_metfe', 'st_metfeerr1', 'st_metfeerr2',
                   'st_teff', 'st_tefferr1', 'st_tefferr2',
                   'st_plx', 'st_plxerr1', 'st_plxerr2',
                   'st_logg', 'st_loggerr1', 'st_loggerr2',
                   'st_mass', 'st_masserr1', 'st_masserr2',
                   'st_spstr']
    fields_nasa = ['pl_name', 'hostname', 'hd_name', 'ra', 'dec', 'rastr', 'decstr',
                   'sy_vmag', 'sy_vmagerr1', 'sy_vmagerr2',
                   'st_met', 'st_meterr1', 'st_meterr2',
                   'st_teff', 'st_tefferr1', 'st_tefferr2',
                   'st_logg', 'st_loggerr1', 'st_loggerr2',
                   'st_mass', 'st_masserr1', 'st_masserr2',
                   'st_spectype', 'pl_refname', 'gaia_id']

    if nasa:
        # Loading NASA exoplanet archive
        exo_all = pd.read_csv('nasaexo.csv',
                              skipinitialspace=True, usecols=fields_nasa)
        # Changing some column names to match exoplanet.EU
        exo_all = exo_all.rename(columns={"pl_name": "name",
                                          "hostname": "star_name",
                                          "sy_vmag": "mag_v",
                                          "sy_vmagerr1": "mag_v_err",
                                          "st_teff": "star_teff",
                                          "st_tefferr1": "star_teff_error_max",
                                          "st_tefferr2": "star_teff_error_min",
                                          "st_met": "star_metallicity",
                                          "st_meterr1": "star_metallicity_error_max",
                                          "st_meterr2": "star_metallicity_error_min"})
    else:
        # Laoding exoplanet.EU
        exo_all = pd.read_csv('exo.csv', skipinitialspace=True, usecols=fields)

    # Remove trailing whitespaces
    exo_all.star_name = exo_all.star_name.str.strip()
    #print(exo_all.star_name)
    output = 'WEBSITE_online_EU-NASA_to_ADD.rdb'
    for i, star in enumerate(stars):
        star = star.strip('\n').replace('K0I', 'KOI')
        print(i, star)
        #exo = exo_all[exo_all.star_name == star]
        exo = exo_all[(exo_all.ra == ra_stars[i]) & (exo_all.dec == dec_stars[i])]
        print(exo)
        exo = exo.reset_index(drop=True)
        #exo = exo.iloc[0]
        next = True
        #print(exo)
        print('')
        print('Star/planet: ' + colored.green(exo['name'][0]))

        try:
            name = exo['star_name'][0]
            if name is None or name == "None":
                name = remove_planet(exo['name'][0])
            print('Star: ' + colored.green(name))
            print("Stars to add: ", len(linesnames))
        except IndexError as e:
            print('')
            puts(colored.red(star) + ' not found. Star added in the file manual.list.')
            print('')
            manual = open('manual.list', "a")
            manual.write(star+'\n')
            manual.close()
            next = False
            # Update the list of new hosts
            with open('names.txt', 'w') as names:
                # if the last star was added so no star is updated
                if i+1==len(stars):
                    names.write('')
                else:
                    for j in stars[i+1:]:
                        names.write(j)
            names.close()
            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

        # if the star is found in the exoplanet.eu
        # or if found in NASA exoplanet archive
        if next:
            print('')
            var = input('Continue? [Y/N]: ')

            if var.upper().strip() == 'Y':

                # Get RA and dec (can be passed for NASA exoplanets)
                ra, dec = ra_stars[i], dec_stars[i]
                c = coord.SkyCoord(ra, dec, unit=(u.degree, u.degree),
                                   frame='icrs')
                RA = list(c.ra.hms)
                RA[0] = str(int(RA[0])).zfill(2)
                RA[1] = str(int(RA[1])).zfill(2)
                RA[2] = str(round(RA[2], 2)).zfill(4)
                if len(RA[2]) == 4:
                    RA[2] += '0'
                RA = "{0} {1} {2}".format(*RA)

                DEC = list(c.dec.dms)
                DEC[0] = str(int(DEC[0])).zfill(2)
                DEC[1] = str(abs(int(DEC[1]))).zfill(2)
                DEC[2] = str(abs(round(DEC[2], 2))).zfill(4)
                if int(DEC[0]) > 0:
                    DEC[0] = '+'+DEC[0]
                if len(DEC[2]) == 4:
                    DEC[2] += '0'
                DEC = "{0} {1} {2}".format(*DEC)

                # Search in Simbad the parallax, Vmag and spectral type
                customSimbad = Simbad()
                # customSimbad.add_votable_fields('plx','plx_error','flux(V)','flux_error(V)','sptype','otype','ids','dist')
                customSimbad.add_votable_fields('plx', 'plx_error', 'flux(V)',
                                                'flux_error(V)', 'sptype',
                                                'otype', 'ids')
                result = customSimbad.query_region(coord.SkyCoord(ra=c.ra,
                                                                  dec=c.dec,
                                                                  frame='icrs'),
                                                   radius='15s')

                empty = 'NULL'
                                
                # Here comes the user interface part...
                puts(colored.black('\nStandard parameters\n'))

                # The metallicity error
                if ~np.isnan(exo.star_metallicity_error_min.values[0]) and ~np.isnan(exo.star_metallicity_error_max.values[0]):
                    errFeH_exo = (abs(exo.star_metallicity_error_min.values[0]) +
                                  abs(exo.star_metallicity_error_max.values[0])) / 2.0
                elif ~np.isnan(exo.star_metallicity_error_min.values[0]):
                    errFeH_exo = abs(exo.star_metallicity_error_min.values[0])
                elif ~np.isnan(exo.star_metallicity_error_max.values[0]):
                    errFeH_exo = abs(exo.star_metallicity_error_max.values[0])
                else:
                    errFeH_exo = np.nan

                # The metallicity
                FeH_exo = exo.star_metallicity.values[0]
                if np.isnan(FeH_exo):
                    puts('The ' + colored.yellow('[Fe/H]'))
                    FeH = variable_assignment(2)
                    puts('The error on ' + colored.yellow('[Fe/H]'))
                    Ferr = variable_assignment(2)
                else:
                    FeH = round(float(FeH_exo), 2)
                    if np.isnan(errFeH_exo):
                        puts('The error on ' + colored.yellow('[Fe/H]'))
                        Ferr = variable_assignment(2)
                    else:
                        Ferr = round(errFeH_exo, 2)
                        
                print('Fe/H: ', FeH, '+-', Ferr)

                # The effective temperature error
                if ~np.isnan(exo.star_teff_error_min.values[0]) and ~np.isnan(exo.star_teff_error_max.values[0]):
                    errTeff_exo = (abs(exo.star_teff_error_min.values[0]) +
                                   abs(exo.star_teff_error_max.values[0])) / 2.0
                elif ~np.isnan(exo.star_teff_error_min.values[0]):
                    errTeff_exo = abs(exo.star_teff_error_min.values[0])
                elif ~np.isnan(exo.star_teff_error_max.values[0]):
                    errTeff_exo = abs(exo.star_teff_error_max.values[0])
                else:
                    errTeff_exo = np.nan

                # The effective temperature
                Teff_exo = exo.star_teff.values[0]
                if np.isnan(Teff_exo):
                    puts('The ' + colored.yellow('Teff'))
                    Teff = variable_assignment(0)
                    puts('The error on ' + colored.yellow('Teff'))
                    Tefferr = variable_assignment(0)
                else:
                    # the Teff is not float
                    Teff = int(Teff_exo)
                    if ~np.isnan(errTeff_exo):
                        Tefferr = int(errTeff_exo)
                    else:
                        puts('The error on ' + colored.yellow('Teff'))
                        Tefferr = variable_assignment(0)

                print('Teff: ', Teff, '+-', Tefferr)

                # NASA database has loggs
                # The logg error and logg
                if nasa:
                    if ~np.isnan(exo.st_loggerr2.values[0]) and ~np.isnan(exo.st_loggerr1.values[0]):
                        errlogg_exo = (abs(exo.st_loggerr2.values[0]) +
                                      abs(exo.st_loggerr1.values[0])) / 2.0
                    elif ~np.isnan(exo.st_loggerr2.values[0]):
                        errlogg_exo = abs(exo.st_loggerr2.values[0])
                    elif ~np.isnan(exo.st_loggerr1.values[0]):
                        errlogg_exo = abs(exo.st_loggerr1.values[0])
                    else:
                        errlogg_exo = np.nan

                    logg_exo = exo.st_logg.values[0]
                    if np.isnan(logg_exo):
                        puts('The ' + colored.yellow('logg'))
                        logg = variable_assignment(2)
                        puts('The error on ' + colored.yellow('logg'))
                        loggerr = variable_assignment(2)
                    else:
                        # logg is a float
                        logg = round(float(logg_exo), 2)
                        if ~np.isnan(errlogg_exo):
                            loggerr = round(errlogg_exo, 2)
                        else:
                            puts('The error on ' + colored.yellow('logg'))
                            loggerr = variable_assignment(0)
                else:
                    # The log g and log g error for EU database
                    puts('The ' + colored.yellow('logg'))
                    logg = variable_assignment(2)
                    puts('The error on ' + colored.yellow('logg'))
                    loggerr = variable_assignment(2)

                print('logg: ', logg, '+-', loggerr)
                #Skipping the mass - will be computed later
                M, Merr = -1,-1

                # The microturbulence number
                puts('The '+colored.yellow('microturbulence'))
                vt = variable_assignment(2)
                puts('The error on '+colored.yellow('microturbulence'))
                vterr = variable_assignment(2)

                if nasa:
                    t = exo.pl_refname[0]
                    print(t)
                    ts = t.replace(" target=ref> ", "00coisinha00").replace('href=',"00coisinha00").replace('target=ref>',"00coisinha00").replace("</a>","").split("00coisinha00")
                    print(ts)
                    print("Link:", ts[1])
                    print("author:", ts[2])
                    author = ts[2].strip()
                    link = ts[1].strip()
                else:
                    # Author and link to ADS
                    puts('Who is the '+colored.yellow('author?'))
                    author = input('> ').strip()
                    if author == '':
                        author = empty
                    puts('Link to article ('+colored.yellow('ADS')+')')
                    link = input('> ').strip()
                    if link == '':
                        link = empty
                # Source flag
                puts(colored.yellow('Source flag'))
                source = input('(0/1) > ')
                if source == '':
                    source = '0'

                V_exo = exo.mag_v.values[0]

                try:
                    # select the star and not the planet,
                    # they have the same coordinates
                    if len(result) > 1:
                        indr = np.where((result['OTYPE'] != 'Planet') &
                                        (result['OTYPE'] != 'Planet?') &
                                        (result['OTYPE'][1] != 'brownD*'))[0][0]
                    else:
                        indr = 0

                    RA = str(result['RA'][indr])[:11]
                    DEC = str(result['DEC'][indr])[:12]

                    # The HD number
                    HD=empty
                    for iname in result['IDS'][indr].split('|'):
                        if iname[:2]=='HD':
                            HD=iname.replace('HD ','')
                    
                    # The V magnitude
                    if type(result['FLUX_V'][indr])!=np.ma.core.MaskedConstant:
                        V=round(float(result['FLUX_V'][indr]), 2)
                        if type(result['FLUX_ERROR_V'][indr])!=np.ma.core.MaskedConstant:    
                            Verr=round(float(result['FLUX_ERROR_V'][indr]), 2)
                        else:
                            print('\nV magnitude = '+str(V))
                            puts('The error on ' + colored.yellow('V magnitude'))
                            Verr = variable_assignment(2)
                            if Verr == '':
                                Verr = 'NULL'                       
                    else:
                        if ~np.isnan(V_exo):
                            V = round(float(V_exo), 2)
                        else:    
                            puts('The ' + colored.yellow('V magnitude'))
                            V = variable_assignment(2)
                            if V == '':
                                V = 'NULL'
                        print('\nV magnitude = '+str(V))
                        puts('The error on ' + colored.yellow('V magnitude'))
                        Verr = variable_assignment(2)
                        if Verr == '':
                            Verr = 'NULL'

                    # The parallax
                    plx, eplx = GAIAplx(RA, DEC)
                    if plx != 'NULL':
                        p = plx
                        perr = eplx
                        pflag = 'GAIADR2'
                    elif type(result['PLX_VALUE'][indr]) != np.ma.core.MaskedConstant:
                        p = round(float(result['PLX_VALUE'][indr]), 2)
                        if type(result['PLX_VALUE'][indr]) != np.ma.core.MaskedConstant:
                            perr = round(float(result['PLX_ERROR'][indr]), 2)
                        else:
                            perr = empty
                        pflag = 'Simbad'
                    else:
                        try:
                            pos = coord.SkyCoord(ra=ra, dec=dec,
                                                 unit=(u.hourangle,u.deg),
                                                 frame='icrs')
                            #AvSF = Schlafly & Finkbeiner 2011 (ApJ 737, 103)
                            tableAv = IrsaDust.get_query_table(pos,
                                                               radius='02d',
                                                               section='ebv',
                                                               timeout=60)
                            Av = tableAv['ext SandF mean'].data[0]
                            Averr = tableAv['ext SandF std'].data[0]
                        except:
                            Av = 0
                            Averr = 0
                        try:    
                            p, perr = [round(x, 2) for x in parallax(Teff,
                                                                     Tefferr,
                                                                     float(logg),
                                                                     float(loggerr),
                                                                     V, Verr, M, Merr, Av, Averr)]
                            pflag = 'Spec'
                        except:
                            p = 'NULL'
                            perr = 'NULL'
                            pflag = 'NULL'

                    # Comments
                    if result['SP_TYPE'][indr] != '' and result['SP_TYPE'][indr][0] == 'M':
                        comment = result['SP_TYPE'][indr]
                    else:
                        puts('Any '+colored.yellow('comments'))
                        puts('E.g. if we have a M dwarf...')
                        comment = input('> ')
                        if comment == '':
                            comment = 'NULL'

                    if nasa:
                        database = "NASA"
                    else:
                        database = "EU"
#                    # Exoplanet database
#                    puts('From which exoplanet database:' + colored.yellow('EU or NASA or EU,NASA'))
#                    database = input('> ')
#                    if database == '':
#                        database = 'NULL'

                except:

                    # The HD number
                    puts('The '+colored.yellow('HD number'))
                    HD = input('> ')
                    if HD == '':
                        HD = 'NULL'

                    # The V magnitude
                    if ~np.isnan(V_exo):
                        V = round(float(V_exo), 2)
                    else:    
                        puts('The ' + colored.yellow('V magnitude'))
                        V = variable_assignment(2)
                    print('\nV magnitude = '+str(V))
                    puts('The error on ' + colored.yellow('V magnitude'))
                    Verr = variable_assignment(2)

                    # The parallax
                    plx, eplx = GAIAplx(RA, DEC)
                    if plx != 'NULL':
                        p = plx
                        perr = eplx
                        pflag = 'GAIADR2'
                    else:
                        try:
                            pos = coord.SkyCoord(ra=RA, dec=DEC,
                                                 unit=(u.hourangle, u.deg),
                                                 frame='icrs')
                            #AvSF = Schlafly & Finkbeiner 2011 (ApJ 737, 103)
                            tableAv = IrsaDust.get_query_table(pos,
                                                               radius='02d',
                                                               section='ebv',
                                                               timeout=60)
                            Av = tableAv['ext SandF mean'].data[0]
                            Averr = tableAv['ext SandF std'].data[0]
                        except:
                            Av = 0
                            Averr = 0
                        try:
                            p, perr = [round(x, 2) for x in parallax(Teff,
                                                                     Tefferr,
                                                                     logg, loggerr,
                                                                     V, Verr, M, Merr, Av, Averr)]
                            pflag = 'Spec'
                            # print p,perr
                        except:
                            p = 'NULL'
                            perr = 'NULL'
                            pflag = 'NULL'

                    # Comments
                    puts('Any '+colored.yellow('comments'))
                    puts('E.g. if we have a M dwarf...')
                    comment = input('> ')
                    if comment == '':
                        comment = 'NULL'
                    if nasa:
                        database = "NASA"
                    else:
                        database = "EU"
#                    # Exoplanet database
#                    puts('From which exoplanet database: ' + colored.yellow('EU or NASA or EU,NASA'))
#                    database = input('> ')
#                    if database == '':
#                        database = 'NULL'

                # Last update
                update = str(time.strftime("%Y-%m-%d"))

################ Add new columns:
                #GAIA DR2

                gaiadr2 = get_gaiadr2(name)
                if gaiadr2 == -1:
                    if nasa:
                        try:
                            gaiadr2 = exo.gaia_id[0].replace("Gaia DR2 ", "")
                        except:
                            gaiadr2 = -1
                    else:    
                        print("GAIA DR2 Not found automatically:")
                        puts('Insert ' + colored.yellow('GAIADR2'))
                        gaiadr2 = input('> ')
                print("GAIADR2:", gaiadr2)

                #Assuming GAIADR3 = GAIADR2
                gaiadr3 = get_gaiadr3(name)
                if gaiadr3 == -1:
                    if nasa:
                        try:
                            gaiadr3 = exo.gaia_id[0].replace("Gaia DR2 ", "")
                        except:
                            gaiadr3 = -1
                    else:
                        print("GAIA EDR3 Not found automatically:")
                        puts('Insert ' + colored.yellow('GAIA EDR3'))
                        gaiadr3 = input('> ')
                print("GAIADR3:", gaiadr3)

                #getting DR3 data:
                if gaiadr2 != "-1":
                    gaia_data, std_g_flux = get_gaiadr3_data(gaiadr2, gaiadr3)
                    if gaia_data['Plx'] == -1:
                        print("Trying to search by name:", name)
                        print("Search by other name?: (Enter if not):")
                        name_in = input('> ')
                        if name_in == "":
                            name_in = name
                        gaia_data, std_g_flux = get_gaiadr3_data(gaiadr2, gaiadr3, name_in)
                    plx_flag = "GAIADR3"
                    try:
                        logg_gaiadr3, erlogg_gaiadr3, masst, emasst, radiust, eradiust = get_logg_dr3_stuff(Teff, Tefferr, logg, loggerr, FeH, Ferr, gaia_data["Gmag"], gaia_data["e_Gmag"], gaia_data["Plx"], gaia_data["e_Plx"])
                    except:
                        logg_gaiadr3, erlogg_gaiadr3, masst, emasst, radiust, eradiust = (-1,-1,-1,-1,-1,-1)
                else:
                    gaia_data  = {'Plx'    : None, 'e_Plx'  : None, 
                                  'FG'     : None, 'e_FG'   : None,  'Gmag'   : None,  'e_Gmag' : None, 
                                  'BPmag'  : None, 'e_BPmag': None,  'RPmag'  : None,  'e_RPmag': None}
                    std_g_flux = None
                    # Manual paralax
                    puts('Manual Paralax ' + colored.yellow('Insert mili arcsec'))
                    gaia_data["Plx"] = read_float()
                    puts('error Paralax ' + colored.yellow('Insert mili arcsec'))
                    gaia_data["e_Plx"] = read_float()
                    puts('source Paralax ' + colored.yellow('DR2, DR3, Hipparcus, Simbad'))
                    plx_flag = input('> ')

                    logg_gaiadr3, erlogg_gaiadr3, masst, emasst, radiust, eradiust  = (None, None, None, None, None, None)
                    if gaia_data["Plx"] > 0:
                        logg_gaiadr3, erlogg_gaiadr3, masst, emasst = loggdr3.logg_mass_hip_iteractive_error(Teff, Tefferr, logg, loggerr, FeH, Ferr, V, Verr, gaia_data["Plx"], gaia_data["e_Plx"])

                print(gaia_data)

                new_row = pd.DataFrame([], columns=columns_SW)
                new_row = new_row.append({"Name": "Empty"}, ignore_index = True)
                new_row["Name"]         = name
                new_row["hd"]           = HD
                new_row["RA"]           = RA
                new_row["DEC"]          = DEC
                new_row["Vmag"]         = V
                new_row["eVmag"]        = Verr
                new_row["Teff"]         = Teff
                new_row["eTeff"]        = Tefferr
                new_row["Logg"]         = logg
                new_row["eLogg"]        = loggerr
                new_row["Vt"]           = vt
                new_row["eVt"]          = vterr
                new_row["[Fe/H]"]       = FeH
                new_row["e[Fe/H]"]      = Ferr
                new_row["Reference"]    = author
                new_row["Link"]         = link
                new_row["SWFlag"]       = source
                new_row["Update"]       = update
                new_row["Comment"]      = comment
                new_row["Database"]     = database
                #...
                new_row["PlxFlag"]      = plx_flag
                new_row["gaia_dr2"]     = str(gaiadr2)
                new_row["gaia_dr3"]     = str(gaiadr3) 

                new_row["Plx"]          = gaia_data["Plx"]
                new_row["ePlx"]         = gaia_data["e_Plx"] 
                new_row["Gmag"]         = gaia_data["Gmag"]
                new_row["eGmag"]        = gaia_data["e_Gmag"]
                new_row["RPmag"]        = gaia_data["RPmag"]
                new_row["eRPmag"]       = gaia_data["e_RPmag"] 
                new_row["BPmag"]        = gaia_data["BPmag"]
                new_row["eBPmag"]       = gaia_data["e_BPmag"]
                new_row["FG"]           = gaia_data["FG"]
                new_row["eFG"]          = gaia_data["e_FG"]
                new_row["G_flux_std_n"] = std_g_flux

                new_row["Logg_gaia"]    = logg_gaiadr3
                new_row["eLogg_gaia"]   = erlogg_gaiadr3
                new_row["Mass_t"]       = masst
                new_row["eMass_t"]      = emasst
                new_row["Radius_t"]     = radiust
                new_row["eRadius_t"]    = eradiust

                new_row["spec_base"]    = None
                new_row["Distance"]     = 1000./new_row["Plx"].astype(float)
                if nasa:
                    new_row["RA_NASA"]      = ra
                    new_row["DEC_NASA"]     = dec
                else:
                    new_row["RA_EU"]        = ra
                    new_row["DEC_EU"]       = dec

                print(new_row)
                new_row.to_csv(output+'.csv', mode='a', header=False, index=False)

                params = [name, HD, RA, DEC, V, Verr, p, perr, pflag,
                          Teff, Tefferr, logg, loggerr, 'NULL', 'NULL',
                          vt, vterr, FeH, Ferr, M, Merr,
                          author, link, source, update, comment, database]
                params = list(map(str, params))

                print('\n'+'\t'.join(params) + '\tNULL')
                # New host information
                with open(output, 'a') as f:
                    f.write('\n'+'\t'.join(params) + '\tNULL')
                f.close()

                # Update the list of new hosts
                linesnames = linesnames[1:]
                if len(linesnames) > 0:
                    filenames = open(filenamesin, "w")
                    for linesname in linesnames:
                        filenames.write(linesname)
                    filenames.close()
                print('')
                print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

            else:
                print('Bye...')
                break
