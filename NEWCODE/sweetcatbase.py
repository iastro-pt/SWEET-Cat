import pandas as pd
import os
from astropy import coordinates as coord
from astropy import units as u
import numpy as np
from astropy.io import votable
import warnings
import urllib.request, urllib.error, urllib.parse
from clint.textui import puts, colored 
import glob
import wget


file_path=os.path.dirname(__file__)

columns_SW = ['Name', 'hd', 'RA', 'DEC', 'Vmag', 'eVmag', 'PlxFlag', 'Teff', 'eTeff',
       'Logg', 'eLogg', 'Vt', 'eVt', '[Fe/H]', 'e[Fe/H]', 'Reference', 'Link',
       'SWFlag', 'Update', 'Comment', 'Database', 'gaia_dr2', 'gaia_dr3',
       'Plx', 'ePlx', 'Gmag', 'eGmag', 'RPmag', 'eRPmag', 'BPmag', 'eBPmag',
       'FG', 'eFG', 'G_flux_std_n', 'Logg_gaia', 'eLogg_gaia', 'Mass_t',
       'eMass_t', 'Radius_t', 'eRadius_t', 'spec_base', 'Distance', 'RA_EU',
       'DEC_EU', 'RA_NASA', 'DEC_NASA', 'Distance_b', 'eDistance_b']

dtype_SW = {'gaia_dr2':'int64','gaia_dr3':'int64'}




def load_SC_database(sweet_database, list_stats=False):
  SC_full_load = pd.read_csv(sweet_database, dtype=dtype_SW)
  if list_stats:
    nph = len(SC_full_load)
    nph_h = len(SC_full[SC_full["SWFlag"]==1])
    print(f"Number of planet hosts: {nph}")
    print(f"Number of homogeneous planet hosts: {nph_h}  ({nph_h/nph*100:.0f}%)")

    print("For brighter planet hosts (Gmag <= 12):")
    SC_bright = SC_full_load[SC_full_load["Gmag"] <= 12]
    nph_b = len(SC_bright)
    nph_b_h = len(SC_bright[SC_bright["SWFlag"]==1])
    print(f"Number of planet hosts: {nph_b}")
    print(f"Number of homogeneous planet hosts: {nph_b_h}  ({nph_b_h/nph_b*100:.0f}%)")

    print("For brighter planet hosts (Gmag <= 12) in the Southern Hemisphere:")
    SC_full_load["DEC_deg"] = [float(d.split()[0])+float(d.split()[1])/60++float(d.split()[2])/3600 for d in SC_full_load["DEC"]] 
    SC_bright_S = SC_full_load[(SC_full_load["Gmag"] <= 12) & (SC_full_load["DEC_deg"] < 0)]
    nph_b = len(SC_bright_S)
    nph_b_h = len(SC_bright_S[SC_bright_S["SWFlag"]==1])
    print(f"Number of planet hosts: {nph_b}")
    print(f"Number of homogeneous planet hosts: {nph_b_h}  ({nph_b_h/nph_b*100:.0f}%)")

    print("For brighter planet hosts (Gmag <= 12) in the Northern Hemisphere:")
    SC_bright_N = SC_full_load[(SC_full_load["Gmag"] <= 12) & (SC_full_load["DEC_deg"] >= 0)]
    nph_b = len(SC_bright_N)
    nph_b_h = len(SC_bright_N[SC_bright_N["SWFlag"]==1])
    print(f"Number of planet hosts: {nph_b}")
    print(f"Number of homogeneous planet hosts: {nph_b_h}  ({nph_b_h/nph_b*100:.0f}%)")

    print("For brighter planet hosts (Gmag <= 12) observed by ESO:")
    SC_bright_E = SC_full_load[(SC_full_load["Gmag"] <= 12) & (SC_full_load["DEC_deg"] <= 30)]
    nph_b = len(SC_bright_E)
    nph_b_h = len(SC_bright_E[SC_bright_E["SWFlag"]==1])
    print(f"Number of planet hosts: {nph_b}")
    print(f"Number of homogeneous planet hosts: {nph_b_h}  ({nph_b_h/nph_b*100:.0f}%)")

#Merge Sweet Cat with exo.eu:


def remove_planet(name):
    """Remove the trailing b, c, d, etc in the stellar name"""
    planets = 'abcdefghijB'
    for planet in planets:
        if name.endswith(' %s' % planet):
            return name[:-2].strip()
        if name.endswith(' A%s' % planet):
            return name[:-1].strip()
        if name.endswith(' B%s' % planet):
            return name[:-1].strip()

    # some exoplanets have .01 or .02 in the name 
    if name.endswith('.01') or name.endswith('.02') or name.endswith('.2'):
        return name[:-3].strip()    
    if  name.endswith(' (bc)'):
        return name[:-4].strip()
    return name

def downloadExoplanet_old():
    """
    Download the table from exoplanetEU and save it to a file (exo.csv).

    Return a pandas DataFrame sorted in 'update'.
    """
    print("Downloading EXO.EU exoplanets...")

    response = urllib.request.urlopen("http://exoplanet.eu/catalog/votable")
    table = response.read()
    with open('exo.xml', 'wb') as f:
        f.write(table)
    """Convert the saved xml file to csv and read with pandas"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        vo = votable.parse('exo.xml', invalid='mask', pedantic=False)
        vo = vo.get_first_table().to_table(use_names_over_ids=True)
        df = vo.to_pandas()

    # Divide the data in Confirmed and not.
    df[df.planet_status == 'Confirmed'].to_csv('exo.csv', index=False)
    print("Downloading EXO.EU exoplanets... DONE!")
    os.system("rm exo.xml")
    return

def downloadExoplanet():
    """
    Download the table from exoplanetEU and save it to a file (exo.csv).

    Return a pandas DataFrame sorted in 'update'.
    """
    print("Downloading EXO.EU exoplanets...")
    filename = "exo.csv"
    if os.path.exists(filename):
      os.remove(filename) # if exist, remove it directly
    urlexo = "https://exoplanet.eu/catalog/csv"
    wget.download(urlexo, out=filename)
    print("Downloading EXO.EU exoplanets... DONE!")
    return



def downloadNasaExoplanetNew():
    """
    Download the table from NASA exoplanet archive
    and save it to a file (nasaexo.csv).
    """
    print("Downloading NASA exoplanets...")
    urlRoot = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
    table = "+from+ps"     # table of planetary systems
    select = "select+"             # which columns to select
    order = "+where+default_flag=1"
    tableFormat = "&format=csv"     # format of the table

    # Choose SweetCat columns
    sweetCatColumns = "*"

    # Full url for the query
    nasa_url = urlRoot  + select + sweetCatColumns + table + order + tableFormat

    # Download the data
    response = urllib.request.urlopen(nasa_url)
    table = response.read()

    # Write the NASA exoplanet archive
    with open('nasaexo.csv', 'wb') as f:
        f.write(table)
    print("Downloading NASA exoplanets... DONE!")
    return

def check_SW_coordinates(SCcsv, EXOcsv, NASAcsv):
  SC = pd.read_csv(SCcsv, dtype=dtype_SW)
  exo = pd.read_csv(EXOcsv)
#  print(np.unique(exo.detection_type))
#  exo = exo[(exo.detection_type == 'Radial Velocity') | (exo.detection_type == 'Radial Velocity, Astrometry') | (exo.detection_type == 'Primary Transit') | (exo.detection_type == 'Astrometry') | (exo.detection_type == 'Default') | (exo.detection_type == 'Other')]
  exo = exo.reset_index()
  nasa = pd.read_csv(NASAcsv)


  coordSC     = coord.SkyCoord(ra=SC['RA'].values,
                               dec=SC['DEC'].values,
                               unit=(u.hourangle, u.deg),
                               frame='icrs')
  coordSCEU   = coord.SkyCoord(ra=SC['RA_EU'].values,
                               dec=SC['DEC_EU'].values,
                               unit=(u.deg, u.deg),
                               frame='icrs')
  coordSCNASA = coord.SkyCoord(ra=SC['RA_NASA'].values,
                               dec=SC['DEC_NASA'].values,
                               unit=(u.deg, u.deg),
                               frame='icrs')
  coordExo    = coord.SkyCoord(ra=exo['ra'].values,
                               dec=exo['dec'].values,
                               unit=(u.deg, u.deg),
                               frame='icrs')
  coordNasa   = coord.SkyCoord(ra=nasa['ra'].values,
                               dec=nasa['dec'].values,
                               unit=(u.deg, u.deg),
                               frame='icrs')

  #check EU
  print("Checking EU:\n")
  for i in range(len(SC)):
    sep = coordSC[i].separation(coordExo).arcsecond
    sepeu = coordSCEU[i].separation(coordExo).arcsecond
    ind = np.where(sep <= np.nanmin(sep))[0]
    indeu = np.where(sepeu <= np.nanmin(sepeu))[0]
    if "EU" in SC.Database[i]:
      if np.isnan(SC.RA_EU[i]) or np.isnan(SC.DEC_EU[i]):
        print("missing eu RA DEC:", i, SC.Name[i])
      else:
        sep_eu = coordSC[i].separation(coordSCEU[i]).arcsecond
        #print(i, SC.Name[i], sep_eu, np.nanmin(sep), sep_eu > np.nanmin(sep))
        #if sep_eu > np.nanmin(sep):
        #  print(i,SC.Name[i], sep_eu, ind[0], exo.name[ind[0]], sepeu[indeu[0]], exo.name[indeu[0]])
        if sepeu[indeu[0]] > 0.01: 
          print("OIOIOI\n",i,SC.Name[i], sep_eu, ind[0], exo.name[ind[0]], sepeu[indeu[0]], exo.name[indeu[0]])

  #check NASA
  print("\n\nChecking NASA:\n")
  for i in range(len(SC)):
    sep = coordSC[i].separation(coordNasa).arcsecond
    sepnasa = coordSCNASA[i].separation(coordNasa).arcsecond
    ind = np.where(sep <= np.min(sep))[0]
    indnasa = np.where(sepnasa <= np.min(sepnasa))[0]
    if "NASA" in SC.Database[i]:
      if np.isnan(SC.RA_NASA[i]) or np.isnan(SC.DEC_NASA[i]):
        print("missing NASA RA DEC:", i, SC.Name[i])
      else:
        sep_nasa = coordSC[i].separation(coordSCNASA[i]).arcsecond
        #print(sep_nasa, np.nanmin(sep))
        #if sep_nasa > np.min(sep):
        #  print(i,SC.Name[i], sep_nasa, nasa.pl_name[ind[0]], sepnasa[indnasa[0]], nasa.pl_name[indnasa[0]])
        if sepnasa[indnasa[0]] > 0.01: 
          print("OIOIOI\n", i,SC.Name[i], sep_nasa, nasa.pl_name[ind[0]], sepnasa[indnasa[0]], nasa.pl_name[indnasa[0]])

  #check GAIA_ID coords:



  return



def askFix_SW_coordinates(SCcsv, EXOcsv, NASAcsv):
  SC = pd.read_csv(SCcsv, dtype=dtype_SW)
  exo = pd.read_csv(EXOcsv)
  exo = exo[(exo.detection_type == 'Radial Velocity') | (exo.detection_type == 'Radial Velocity, Astrometry') | (exo.detection_type == 'Primary Transit') | (exo.detection_type == 'Astrometry') | (exo.detection_type == 'Default') | (exo.detection_type == 'Other')]
  exo = exo.reset_index()
  nasa = pd.read_csv(NASAcsv)


  coordSC     = coord.SkyCoord(ra=SC['RA'].values,
                               dec=SC['DEC'].values,
                               unit=(u.hourangle, u.deg),
                               frame='icrs')
  coordSCEU   = coord.SkyCoord(ra=SC['RA_EU'].values,
                               dec=SC['DEC_EU'].values,
                               unit=(u.deg, u.deg),
                               frame='icrs')
  coordSCNASA = coord.SkyCoord(ra=SC['RA_NASA'].values,
                               dec=SC['DEC_NASA'].values,
                               unit=(u.deg, u.deg),
                               frame='icrs')
  coordExo    = coord.SkyCoord(ra=exo['ra'].values,
                               dec=exo['dec'].values,
                               unit=(u.deg, u.deg),
                               frame='icrs')
  coordNasa   = coord.SkyCoord(ra=nasa['ra'].values,
                               dec=nasa['dec'].values,
                               unit=(u.deg, u.deg),
                               frame='icrs')

  #check EU
  print("Checking EU:\n")
  for i in range(len(SC)):
    sep = coordSC[i].separation(coordExo).arcsecond
    sepeu = coordSCEU[i].separation(coordExo).arcsecond
    ind = np.where(sep <= np.nanmin(sep))[0]
    indeu = np.where(sepeu <= np.nanmin(sepeu))[0]
    if "EU" in SC.Database[i]:
      if np.isnan(SC.RA_EU[i]) or np.isnan(SC.DEC_EU[i]):
        print("missing eu RA DEC:", i, SC.Name[i])
      else:
        sep_eu = coordSC[i].separation(coordSCEU[i]).arcsecond
        #print(i, SC.Name[i], sep_eu, np.nanmin(sep), sep_eu > np.nanmin(sep))
        #if sep_eu > np.nanmin(sep):
        #  print(i,SC.Name[i], sep_eu, ind[0], exo.name[ind[0]], sepeu[indeu[0]], exo.name[indeu[0]])
        if sepeu[indeu[0]] > 0.01: 
          print("OIOIOI\n",i,SC.Name[i], sep_eu, ind[0], exo.name[ind[0]], sepeu[indeu[0]], exo.name[indeu[0]])
          var = input('Update Coordinates? [Y/N]: ')
          if var.upper() == 'Y':
            print("Need to update:",exo.name[ind[0]],exo.ra[ind[0]], exo.dec[ind[0]])
            print("BEFORE:", SC.RA_EU[i],SC.DEC_EU[i])
            SC.RA_EU[i] = exo.ra[ind[0]]
            SC.DEC_EU[i] = exo.dec[ind[0]]
            print("AFTER", SC.RA_EU[i],SC.DEC_EU[i])

  SC.to_csv('SWEETCAT_Dataframe_updatedCoords.csv', index=False)


def check_missing_SweetCat_ExoEU2(SCcsv, EXOcsv):
  exo = pd.read_csv(EXOcsv)
  SC = pd.read_csv(SCcsv, dtype=dtype_SW)
  print(np.unique(exo.detection_type))
  allow_detection_list = ['Radial Velocity', 'Primary Transit', 'Astrometry', 'Imaging, Astrometry', 'Imaging, Other, Astrometry',
 'Imaging, Primary Transit', 'Imaging, Radial Velocity, Astrometry', 'Other, Radial Velocity', 'Primary Transit, Astrometry', 
 'Primary Transit, Radial Velocity', 'Primary Transit, TTV', 'Radial Velocity, Astrometry',  'Radial Velocity, Imaging',
 'Radial Velocity, Other', 'Radial Velocity, Primary Transit', 'Radial Velocity, Timing', 'Timing, Astrometry']
  exo = exo[(exo.detection_type.apply(lambda x: x in allow_detection_list))]
  print(np.unique(exo.detection_type))
  exo =  exo[(exo.planet_status == 'Confirmed')  & (exo.star_name.notna())]
#  exo = exo[(exo.detection_type == 'Radial Velocity') | (exo.detection_type == 'Primary Transit') | (exo.detection_type == 'Astrometry')]
  exo = exo.reset_index()
  coordExo = coord.SkyCoord(ra=exo['ra'].values,
                            dec=exo['dec'].values,
                            unit=(u.deg, u.deg),
                            frame='icrs')
  coordSC = coord.SkyCoord(ra=SC['RA'].values,
                           dec=SC['DEC'].values,
                           unit=(u.hourangle, u.deg),
                           frame='icrs')
  coordSCEU = coord.SkyCoord(ra=SC['RA_EU'].values,
                            dec=SC['DEC_EU'].values,
                            unit=(u.deg, u.deg),
                            frame='icrs')
#  for p in range(1000):
  names_eu_to_add = []
  exo_iloc = []
  for p in range(len(exo)):
    if np.isnan(coordExo[p].ra) or np.isnan(coordExo[p].dec):
      print(exo.name[p], "has problems with star's coordinate")
      continue
    sep = coordExo[p].separation(coordSCEU).arcsecond
    sep = np.array([1e10 if s is None else s for s in sep])
    ind = np.where(sep <= np.nanmin(sep))[0]
    if sep[ind[0]] > 0.01:
      print(p, exo.name[p],ind[0], SC.Name[ind[0]], ind, sep[ind])
      names_eu_to_add.append(exo.name[p])
      exo_iloc.append(p)
  names_eu_to_add = [remove_planet(n) for n in names_eu_to_add]
  names_eu_to_add_u = np.unique(names_eu_to_add)
  fileadd = open("names.txt", 'w')
  for n in names_eu_to_add_u:
    p = exo_iloc[names_eu_to_add.index(n)]
    ####ATENTION Including a limit of maximum mass from exoplanet.eu
    if exo['mass'][p] < 20 or exo['mass'][p] is None:
      fileadd.write(str(n)+"\t"+str(exo['ra'][p])+"\t"+str(exo['dec'][p])+"\t"+str(exo['star_teff'][p])+"\t"+str(exo['mag_v'][p])+"\t"+str(exo['mass'][p])+"\n")

  fileadd.close()
  print(names_eu_to_add)


def check_missing_SweetCat_NASA(SCcsv, NASAcsv):
  nasa = pd.read_csv(NASAcsv)
  nasa = nasa[(nasa.discoverymethod == 'Radial Velocity') | (nasa.discoverymethod == 'Transit') | (nasa.discoverymethod == 'Astrometry')]
  nasa = nasa.reset_index()
  SC = pd.read_csv(SCcsv, dtype=dtype_SW)

  coordNasa   = coord.SkyCoord(ra=nasa['ra'].values,
                               dec=nasa['dec'].values,
                               unit=(u.deg, u.deg),
                               frame='icrs')
  coordSC     = coord.SkyCoord(ra=SC['RA'].values,
                               dec=SC['DEC'].values,
                               unit=(u.hourangle, u.deg),
                               frame='icrs')
  coordSCNASA = coord.SkyCoord(ra=SC['RA_NASA'].values,
                               dec=SC['DEC_NASA'].values,
                               unit=(u.deg, u.deg),
                               frame='icrs')
#  for p in range(1000):
  names_nasa_to_add = []
  nasa_iloc = []
  for n in range(len(nasa)):
    sep = coordNasa[n].separation(coordSCNASA).arcsecond
    sep = np.array([1e10 if s is None else s for s in sep])
    ind = np.where(sep <= np.nanmin(sep))[0]
    if sep[ind[0]] > 0.01:
      print(n, nasa.pl_name[n], nasa.hostname[n], ind[0], SC.Name[ind[0]], ind, sep[ind])
      names_nasa_to_add.append(nasa.hostname[n])
      nasa_iloc.append(n)

  names_nasa_to_add_u = np.unique(names_nasa_to_add)
  fileadd = open("namesnasa.txt", 'w')
  for n in names_nasa_to_add_u:
    p = nasa_iloc[names_nasa_to_add.index(n)]
    fileadd.write(str(n)+"\t"+str(nasa['ra'][p])+"\t"+str(nasa['dec'][p])+"\n")

  fileadd.close()
  print(names_nasa_to_add_u)



def add_nasa_flag_SC(fileSW):
  SC = pd.read_csv(fileSW, dtype=dtype_SW)
  nasanames = pd.read_csv("namesnasa.txt", delimiter="\t", names=["star","ra","dec"])
  coordsnasanames = coord.SkyCoord(ra=nasanames['ra'].values,
                               dec=nasanames['dec'].values,
                               unit=(u.deg, u.deg), frame='icrs')
  coordSC     = coord.SkyCoord(ra=SC['RA'].values,
                               dec=SC['DEC'].values,
                               unit=(u.hourangle, u.deg),
                               frame='icrs')
  droprows = []
  for n in range(len(nasanames)):
    sep = coordsnasanames[n].separation(coordSC).arcsecond
    ind = np.where(sep <= np.nanmin(sep))[0]
    print(ind)
    print(n,len(nasanames), nasanames["star"][n],SC.Name[ind[0]], sep[ind[0]], SC.Database[ind[0]], ind[0], SC.RA_NASA[ind[0]])
    if sep[ind[0]] < 10:
      val = input(colored.green("Accept? ([Y]/N): "))
      if val == "" or val.upper() == "Y":
        val = True
      else:
        val = False
    else:
      val = input(colored.red("Accept? (Y/[N]): "))
      if val == "" or val.upper() == "N":
        val = False
      else:
        val = True
    print(val)
    if val:
      SC.at[ind[0],"Database"] = "EU,NASA"
      droprows.append(n)
      SC.at[ind[0],"RA_NASA"] = coordsnasanames[n].ra.value
      SC.at[ind[0],"DEC_NASA"] = coordsnasanames[n].dec.value
    print(nasanames["star"][n],SC.Name[ind[0]], sep[ind[0]], SC.Database[ind[0]])
    print(len(droprows))
  nasanames=nasanames.drop(droprows)
  SC.to_csv("update.csv", index=False)
  nasanames.to_csv("namesnasa_c.txt", index=False, sep="\t", header=False)


def add_EU_flag_SC(fileSW):
  SC = pd.read_csv(fileSW, dtype=dtype_SW)
  eunames = pd.read_csv("names.txt", delimiter="\t", names=["star","ra","dec", "teff", "magv", "mass"])
  coordseunames = coord.SkyCoord(ra=eunames['ra'].values,
                               dec=eunames['dec'].values,
                               unit=(u.deg, u.deg), frame='icrs')
  coordSC     = coord.SkyCoord(ra=SC['RA'].values,
                               dec=SC['DEC'].values,
                               unit=(u.hourangle, u.deg),
                               frame='icrs')
  droprows = []
  for n in range(len(eunames)):
    sep = coordseunames[n].separation(coordSC).arcsecond
    ind = np.where(sep <= np.nanmin(sep))[0]
    print(ind)
    print(n,len(eunames), eunames["star"][n], SC.Name[ind[0]], sep[ind[0]], SC.Database[ind[0]], ind[0], SC.RA_EU[ind[0]])
    if sep[ind[0]] < 10:
      val = input(colored.green("Accept? ([Y]/N): "))
      if val == "" or val.upper() == "Y":
        val = True
      else:
        val = False
    else:
      val = input(colored.red("Accept? (Y/[N]): "))
      if val == "" or val.upper() == "N":
        val = False
      else:
        val = True
    print(val)
    if val:
      SC.at[ind[0],"Database"] = "EU,NASA"
      droprows.append(n)
      SC.at[ind[0],"RA_EU"] = coordseunames[n].ra.value
      SC.at[ind[0],"DEC_EU"] = coordseunames[n].dec.value
    print(eunames["star"][n],SC.Name[ind[0]], sep[ind[0]], SC.Database[ind[0]])
    print(len(droprows))
  eunames=eunames.drop(droprows)
  SC.to_csv("update.csv", index=False)
  eunames.to_csv("names_c.txt", index=False, sep="\t", header=False)




def check_spectra_present(fileSW, base_dir = "spectra_database/"):
  subfolders = ["combined_spec/", "old_spec/", "added_spectra/", "added_spectra/mstar_spectra/SWEETCAT_PARAMETERS_SPECTRA/"]
  SC = pd.read_csv(fileSW, dtype=dtype_SW)
  SC_miss_specBase = SC[(SC.SWFlag == 1) & (SC.spec_base.isnull())]
  if len(SC_miss_specBase) > 0:
    print("List of homogeneous without spec_base")
    print(SC_miss_specBase[["Name", "spec_base", "SWFlag"]])
  else:
    print("Table with specbase in place")

  ## Checking spec fits in drive
  spec_bases = np.transpose(SC[SC.SWFlag == 1][["spec_base"]].values)[0]
  for i in range(len(spec_bases)):
    fits_file = spec_bases[i]+".fits"
    allgood = True
    files = []
    for sfolder in subfolders:
      files = files + glob.glob(base_dir+sfolder+fits_file)
    if len(files) < 1:
      print("problem: %d - %s  not found" % (i, fits_file))
      allgood = False
  if allgood:
    print("All spectra in Database")


### Main program:
def main():

  fileSW = 'webpage_html/download/SWEETCAT_Dataframe.csv'

## Check spectra in database:
#  check_spectra_present(fileSW)
#  return

##0 Download exo and nasaexo:
  download_planets = False

  if download_planets:
    downloadExoplanet()
    downloadNasaExoplanetNew()


##1: Check coordinates:
  #check_SW_coordinates(fileSW, "exo.csv", "nasaexo.csv")
  #askFix_SW_coordinates(fileSW, "exo.csv", "nasaexo.csv")
  #return

##2: Check EU missing:
  check_missing_SweetCat_ExoEU2(fileSW, "exo.csv")
  #return

#ADD missing and Update file

##3: Check NASA missing:
  #check_missing_SweetCat_NASA(fileSW, "nasaexo.csv")
  #return

##4: Add missing NASA/EU flag (namesnasa.txt):
  #add_nasa_flag_SC(fileSW)
  #add_EU_flag_SC(fileSW)
  #return

## use the namesnasa_c and update "update.csv"

## replace update with the dataframe in the end









if __name__ == "__main__":
    main()