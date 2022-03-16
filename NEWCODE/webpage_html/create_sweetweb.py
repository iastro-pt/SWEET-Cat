import pandas as pd
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import BoxStyle
import astropy.coordinates as coord
from astropy import units as u
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter
from datetime import date

px  = 1/plt.rcParams['figure.dpi']  # pixel in inches


columns_SW = ['Name', 'hd', 'RA', 'DEC', 'Vmag', 'eVmag', 'PlxFlag', 'Teff', 'eTeff',
       'Logg', 'eLogg', 'Vt', 'eVt', '[Fe/H]', 'e[Fe/H]', 'Reference', 'Link',
       'SWFlag', 'Update', 'Comment', 'Database', 'gaia_dr2', 'gaia_dr3',
       'Plx', 'ePlx', 'Gmag', 'eGmag', 'RPmag', 'eRPmag', 'BPmag', 'eBPmag',
       'FG', 'eFG', 'G_flux_std_n', 'Logg_gaia', 'eLogg_gaia', 'Mass_t',
       'eMass_t', 'Radius_t', 'eRadius_t', 'spec_base', 'Distance', 'RA_EU',
       'DEC_EU', 'RA_NASA', 'DEC_NASA', 'Distance_b', 'eDistance_b']

dtype_SW = {'gaia_dr2':'int64','gaia_dr3':'int64'}

dic = {'name':'pl_name','planet_status':'','mass':'pl_bmassj',
       'mass_error_min':'pl_bmassjerr2','mass_error_max':'pl_bmassjerr1',
       'mass_sini':'pl_msinij','mass_sini_error_min':'pl_msinijerr2',
       'mass_sini_error_max':'pl_msinijerr1',
       'radius':'pl_radj', 'radius_error_min':'pl_radjerr2', 
       'radius_error_max':'pl_radjerr1', 'orbital_period':'pl_orbper', 
       'orbital_period_error_min':'pl_orbpererr2', 'orbital_period_error_max':'pl_orbpererr1', 
       'semi_major_axis':'pl_orbsmax', 'semi_major_axis_error_min':'pl_orbsmaxerr2', 
       'semi_major_axis_error_max':'pl_orbsmaxerr1', 'eccentricity':'pl_orbeccen', 
       'eccentricity_error_min':'pl_orbeccenerr2', 'eccentricity_error_max':'pl_orbeccenerr1', 
       'inclination':'pl_orbincl', 'inclination_error_min':'pl_orbinclerr2', 
       'inclination_error_max':'pl_orbinclerr1', 'angular_distance':'', 'discovered':'', 
       'updated':'rowupdate', 'omega':'', 'omega_error_min':'', 
       'omega_error_max':'', 'tperi':'', 'tperi_error_min':'', 'tperi_error_max':'', 
       'tconj':'', 'tconj_error_min':'', 'tconj_error_max':'', 'tzero_tr':'', 'tzero_tr_error_min':'',
       'tzero_tr_error_max':'', 'tzero_tr_sec':'', 'tzero_tr_sec_error_min':'',
       'tzero_tr_sec_error_max':'', 'lambda_angle':'', 'lambda_angle_error_min':'',
       'lambda_angle_error_max':'', 'impact_parameter':'',
       'impact_parameter_error_min':'', 'impact_parameter_error_max':'', 'tzero_vr':'',
       'tzero_vr_error_min':'', 'tzero_vr_error_max':'', 'k':'', 'k_error_min':'',
       'k_error_max':'', 'temp_calculated':'', 'temp_calculated_error_min':'',
       'temp_calculated_error_max':'', 'temp_measured':'', 'hot_point_lon':'',
       'geometric_albedo':'', 'geometric_albedo_error_min':'',
       'geometric_albedo_error_max':'', 'log_g':'', 'publication':'', 'detection_type':'discoverymethod',
       'mass_detection_type':'', 'radius_detection_type':'', 'alternate_names':'',
       'molecules':'', 'star_name':'hostname', 'ra':'ra', 'dec':'dec', 'mag_v':'', 'mag_i':'', 'mag_j':'',
       'mag_h':'', 'mag_k':'', 'star_distance':'sy_dist', 'star_distance_error_min':'sy_disterr2',
       'star_distance_error_max':'sy_disterr1', 'star_metallicity':'',
       'star_metallicity_error_min':'', 'star_metallicity_error_max':'', 'star_mass':'st_mass',
       'star_mass_error_min':'st_masserr2', 'star_mass_error_max':'st_masserr1', 'star_radius':'st_rad',
       'star_radius_error_min':'st_raderr2', 'star_radius_error_max':'st_raderr1', 'star_sp_type':'',
       'star_age':'', 'star_age_error_min':'', 'star_age_error_max':'', 'star_teff':'st_teff',
       'star_teff_error_min':'st_tefferr2', 'star_teff_error_max':'st_tefferr1', 'star_detected_disc':'',
       'star_magnetic_field':'', 'star_alternate_names':''}



#Create a Web-Page from the SC_full DataFrame
def create_html_table(SC_full, filename_web_out = "index.html"):

  print(SC_full)
  print(SC_full.columns)

  columns_web        = ["Name", "RA"  , "DEC" , "gaia_dr3", "Gmag"   , "Plx"    , "Distance"   , "Teff", "eTeff", "Logg", "eLogg", "[Fe/H]", "e[Fe/H]", "Vt"  , "eVt" , "Logg_gaia", "eLogg_gaia", "Mass_t", "eMass_t", "Radius_t", "eRadius_t", "SWFlag", "Reference"]
  columns_web_format = ["%s"  , "%s"  , "%s"  , "%s"      , "%5.2f"  , "%5.2f"  , "%6.1f"      , "%d"  , "%d"   , "%s"  , "%s"   , "%s"    , "%s"     , "%s"  , "%s"  , "%5.3f"    , "%5.3f"      , "%5.3f" , "%5.3f"   , "%7.3f"   , "%7.3f"     , "%s"    , "%s"       ]

  columns_number = ["Gmag", "Plx" , "Distance", "Teff", "eTeff", "Logg", "eLogg", "e[Fe/H]", "Vt"  , "eVt" , "logg_gaia", "e_logg_gaia", "mass_t", "e_mass_t", "radius_t", "e_radius_t"]

  SC_html = SC_full[columns_web]

  filename_web_base = "BaseWeb.html"
  fileweb = open(filename_web_base, "r")
  strlines = fileweb.readlines()
  fileweb.close()

  iline_rows =  strlines.index("<!-- INSERT ROWS HERE !-->\n")

  strlines_out = strlines[:iline_rows]


  for i in range(len(SC_html)):
#  for i in range(100):
    link_simbad = "http://simbad.u-strasbg.fr/simbad/sim-id?Ident=GAIA+DR2+"+str(SC_full["gaia_dr2"][i])+"&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id"
    strline = "<tr> "
    for ir, r in enumerate(columns_web):
      if r == "Name":
        strline += '<td> <a href="'+ link_simbad +'" target="_blank">' + columns_web_format[ir] % (SC_html[r][i]) + '</a></td>' 
      elif r == "Reference":
        strline += '<td> <a href="'+ str(SC_full["Link"][i]) +'" target="_blank">' + columns_web_format[ir] % (SC_html[r][i]) + '</a></td>' 
      else:
        if r in columns_number and float(SC_html[r][i]) < 0:
          strline += '<td> - </td>'
        else:
          try:
            if columns_web_format[ir] == "%s":
              strline += '<td> '+ columns_web_format[ir] % (SC_html[r][i]) + "</td>"
            else:
              strline += '<td> '+ columns_web_format[ir] % (float(SC_html[r][i])) + "</td>"
          except:
              strline += '<td> - </td>'

    strline += "</tr>\n"
    print(strline)
    strline = strline.replace("nan", "-")
    strlines_out.append(strline)
  strlines_out += strlines[iline_rows+1:]
  
  filewebout = open(filename_web_out, "w")
  for line in strlines_out:
    filewebout.write(line)
  filewebout.close()


def SC_stats(SC_full_load, plot_stats=False):
  nph = len(SC_full_load)
  nph_h = len(SC_full_load[SC_full_load["SWFlag"]==1])
  print(f"Number of planet hosts: {nph}")
  print(f"Number of homogeneous planet hosts: {nph_h}  ({nph_h/nph*100:.0f}%)")

  print("For brighter planet hosts (Gmag <= 12):")
  SC_bright = SC_full_load[SC_full_load["Gmag"] <= 12]
  nph_b_1 = len(SC_bright)
  nph_b_h_1 = len(SC_bright[SC_bright["SWFlag"]==1])
  print(f"Number of planet hosts: {nph_b_1}")
  print(f"Number of homogeneous planet hosts: {nph_b_h_1}  ({nph_b_h_1/nph_b_1*100:.0f}%)")

  print("For brighter planet hosts (Gmag <= 12) in the Southern Hemisphere:")
  SC_full_load["DEC_deg"] = [float(d.split()[0])+float(d.split()[1])/60++float(d.split()[2])/3600 for d in SC_full_load["DEC"]] 
  SC_bright_S = SC_full_load[(SC_full_load["Gmag"] <= 12) & (SC_full_load["DEC_deg"] < 0)]
  nph_b_s = len(SC_bright_S)
  nph_b_h_s = len(SC_bright_S[SC_bright_S["SWFlag"]==1])
  print(f"Number of planet hosts: {nph_b_s}")
  print(f"Number of homogeneous planet hosts: {nph_b_h_s}  ({nph_b_h_s/nph_b_s*100:.0f}%)")

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

  if plot_stats:
    cm = plt.cm.get_cmap('YlOrBr')
    fig = plt.figure(figsize=(800*px, 250*px))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    [ax.axis('equal') for ax in [ax1,ax2,ax3]]

#    color_homo = "lightgreen"
#    color_lite = "#FFA500"
    color_homo = "deepskyblue"
    color_lite = "darkorange"

    colors = [color_homo, color_lite]

    #explosion
    explode = (0.05,0.05)

    def make_autopct(values):
      def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.1f}%\n({v:d})'.format(p=pct,v=val)
      return my_autopct

    props = dict(boxstyle='round', facecolor='cyan', alpha=0.5)

    #ax1.set_title("Full Table")
    labels = ['Homogeneous', 'Literature']
    values = [nph_h, nph-nph_h]
    ax1.pie(values,colors=colors, autopct=make_autopct(values), wedgeprops=dict(width=0.6), startangle=90, explode=explode, textprops={'fontsize': 8},pctdistance=0.7)
    ax1.text(0.5, 0.99, "Full Table", transform=ax1.transAxes, fontsize=10, verticalalignment='bottom', horizontalalignment='center', bbox=props)
#    centre_circle1 = plt.Circle((0,0),0.70,fc='cyan')
#    ax1.add_artist(centre_circle1)

    #ax2.set_title("Bright hosts (G<12)")
    labels = ['Homogeneous', 'Literature']
    values = [nph_b_h_1, nph_b_1-nph_b_h_1]
    ax2.pie(values,colors=colors, autopct=make_autopct(values), wedgeprops=dict(width=0.6), startangle=0, explode=explode, textprops={'fontsize': 8},pctdistance=0.7)
    ax2.text(0.5, 0.99, "Bright stars (G < 12)", transform=ax2.transAxes, fontsize=10, verticalalignment='bottom', horizontalalignment='center', bbox=props)
#    centre_circle2 = plt.Circle((0,0),0.70,fc='cyan')
#    ax2.add_artist(centre_circle2)

    #ax3.set_title("Southern Bright hosts")
    labels = ['Homogeneous', 'Literature']
    values = [nph_b_h_s, nph_b_s-nph_b_h_s]
    ax3.pie(values,colors=colors, autopct=make_autopct(values), wedgeprops=dict(width=0.6), startangle=0, explode=explode, textprops={'fontsize': 8},pctdistance=0.7)
    ax3.text(0.5, 0.99, "Southern Bright Stars (Dec < 0)", transform=ax3.transAxes, fontsize=10, verticalalignment='bottom', horizontalalignment='center', bbox=props)
#    centre_circle3 = plt.Circle((0,0),0.70,fc='cyan')
#    ax3.add_artist(centre_circle3)
    props_label = dict(facecolor=color_homo, edgecolor=color_homo)
    props_label2 = dict(facecolor=color_lite, edgecolor=color_lite)
    ax2.text(0.45,0.01, "Homogeneous (SWFlag = 1)", fontsize=8, transform=ax2.transAxes, bbox=props_label, verticalalignment='top', horizontalalignment='right')
    ax2.text(0.55,0.01, "Literature (SWFlag = 0)", fontsize=8, transform=ax2.transAxes, bbox=props_label2, verticalalignment='top', horizontalalignment='left')

    fig.savefig("stat_chart.png")

    plt.show()


def plot_sky_SC_chart(SC):
  SCh = SC[SC.SWFlag == 1]

  coordSC     = coord.SkyCoord(ra=SC['RA'].values,
                dec=SC['DEC'].values,
                unit=(u.hourangle, u.deg),
                frame='icrs')
  coordSCh     = coord.SkyCoord(ra=SCh['RA'].values,
                dec=SCh['DEC'].values,
                unit=(u.hourangle, u.deg),
                frame='icrs')

  ra = coordSC.ra
  dec = coordSC.dec
  ra = ra.wrap_at(180*u.degree)
  rah = coordSCh.ra
  dech = coordSCh.dec
  rah = rah.wrap_at(180*u.degree)

  fig = plt.figure(figsize=(12,8))
  fig.patch.set_facecolor('cyan')
  ax = fig.add_subplot(111, projection="mollweide")
#  ax.scatter(ra.radian, dec.radian, c="darkorange",s=40, alpha=0.75, cmap=cm.brg, label="SWFlag = 0")
#  ax.scatter(rah.radian, dech.radian, c="deepskyblue",s=30, alpha=0.75, cmap=cm.brg, label="SWFlag = 1")
  ax.scatter(ra.radian, dec.radian, c="darkorange",s=40, alpha=0.75, cmap=cm.brg, label="Literature Values")
  ax.scatter(rah.radian, dech.radian, c="deepskyblue",s=30, alpha=0.75, cmap=cm.brg, label="Homogeneous Values")
  ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
  ax.set_facecolor('#0b1925')
  ax.grid(color='white', linestyle='--', linewidth=0.5)
  [t.set_color('white') for t in ax.xaxis.get_ticklines()]
  [t.set_color('white') for t in ax.xaxis.get_ticklabels()]
  [t.set_fontweight('bold') for t in ax.xaxis.get_ticklabels()]
  [t.set_color('k') for t in ax.yaxis.get_ticklines()]
  [t.set_color('k') for t in ax.yaxis.get_ticklabels()]
  [t.set_fontweight('bold') for t in ax.yaxis.get_ticklabels()]
  ax.grid(True)
  ax.legend(loc='upper right')
  today = date.today()
  #Textual month, day and year  
  d2 = today.strftime("%d/%m/%Y")
  ax.text(0.78, 0.10, "Generated on: %s" % (d2), transform=fig.transFigure);
  fig.tight_layout()
  fig.savefig("SC_sky_chart.png")
  plt.show()





def plot_feh_mass_period_all(datafile, homog=1, field_mass='mass_taken', field_period='orbital_period'):
    cm = plt.cm.get_cmap('YlOrBr')
    color_b = (205./255.,205./255., 0./255.)
    color_r = (170./255., 110./255., 40./255.)

    data_rdb = pd.read_csv(datafile)
    data_homog = data_rdb[data_rdb['swsource'] == 1]
    print (len(data_homog), len(data_rdb))

    if homog == 1:
        data_set = data_homog
        print ("homogeneous set")
    else:
        data_set = data_rdb
        print ("full set")
    mass_earth = 1. / 317.8  #jupiter masses
    mass_split = 30. * mass_earth

#    data_set = data_set[np.where( (data_set['teff'] > 4500) & (~np.isnan(data_set['feh'])) & (~np.isnan(data_set[field_mass])) & (~np.isnan(data_set[field_period])) )]
#    data_set = data_set[np.where( (data_set['teff'] > 4500) & (~np.isnan(data_set['feh'])) & (~np.isnan(data_set[field_mass])) & (data_set['detection_type'] == "Radial Velocity") )]
#    data_set = data_set[np.where( (data_set['teff'] > 4500) & (~np.isnan(data_set['feh'])) & (~np.isnan(data_set[field_mass])) & (data_set['detection_type'] == "Primary Transit") )]

    data_set = data_set[(data_set['swteff'] > 4500) & (~np.isnan(data_set['swfeh'])) & (~np.isnan(data_set[field_mass])) & (~np.isnan(data_set[field_period]))]

    print(len(data_set))

    data_lowm = data_set[data_set[field_mass] < mass_split]
    data_highm = data_set[data_set[field_mass] >= mass_split]

    fig = plt.figure(1, figsize=(12, 8))
    ax1 = plt.axes([0.1, 0.1, 0.60, 0.85])
    ax2 = plt.axes([0.72, 0.1, 0.23, 0.85])

    sc = ax1.scatter(data_set['swfeh'], data_set[field_period].astype(float), c=data_set[field_mass].astype(float)/mass_earth,cmap=cm, s=50, ec='k', norm= LogNorm(np.nanmin(data_set[field_mass])/mass_earth, np.nanmax(data_set[field_mass])/mass_earth))
    nullfmt = NullFormatter()         # no labels
    ax1.set_yscale('log')
    cbar = plt.colorbar(sc)
    ax1.set_xlabel('Stellar metallicity - [Fe/H] - (dex)')
    ax1.set_ylabel('Period (days)')
    bins= np.logspace(-1,5,45)
    ax2.set_yscale('log')




    ax2.hist(data_highm[field_period].astype(float), bins=bins, orientation='horizontal', color = color_r, alpha=0.5,ec='black', label='HMP (n: %d) ' % len(data_highm[field_period]))
    ax2.hist(data_lowm[field_period].astype(float), bins=bins, orientation='horizontal', color = color_b, alpha=0.5,ec='black', label='LMP: (n: %d) ' % len (data_lowm[field_period]))
    ax2.set_xlabel('Number of planets')
    ax2.yaxis.set_major_formatter(nullfmt)
    ax2.legend()
#    cbar.set_label(r'Planet Mass (M$_\mathrm{Jupiter}$)', rotation=90)
    cbar.set_label('Planet Mass ($M_{\oplus}$)', rotation=90)
    
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
        ax.set_ylim([10**-1,10**5])
    today = date.today()
    #Textual month, day and year  
    d2 = today.strftime("%d/%m/%Y")
    ax.text(0.15, 0.90, "Generated on: %s" % (d2), transform=fig.transFigure);
    plt.show()
    fig.savefig("mass_metallicity_plane_all.png")

def plot_histograms_parameters(SC):
  fig = plt.figure(figsize=(12,10))
  ax1 = fig.add_subplot(311)
  ax2 = fig.add_subplot(312)
  ax3 = fig.add_subplot(313)

  ax1.hist(SC.Teff, bins = np.arange(2500,7500,100), alpha=0.5,ec='black', color ="deepskyblue")
  ax2.hist(SC["[Fe/H]"], bins = np.arange(-1,0.6,0.05), alpha=0.5,ec='black', color ="deepskyblue")
  ax3.hist(SC["Mass_t"], bins = np.arange(0,3,0.075), alpha=0.5,ec='black', color ="deepskyblue")

  xtitles = ["Teff [K]", "[Fe/H] [dex]", "Stellar Mass [M_sun]"]
  for i, ax in enumerate([ax1, ax2, ax3]):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
    ax.set_ylabel("Number of stars")
    ax.set_xlabel(xtitles[i])
  today = date.today()
  #Textual month, day and year  
  d2 = today.strftime("%d/%m/%Y")
  ax.text(0.78, 0.04, "Generated on: %s" % (d2), transform=fig.transFigure);
  fig.tight_layout()
  plt.show()
  fig.savefig("SC_distributions.png")


#Merge Sweet Cat with exo.eu:
###################################################################################################333



def get_planet_propreties(load = True):


  def get_planets_for_SC_star(SCind):
    outres = []
    if "EU" in SC.Database[SCind]:
      sep = coordSCEU[SCind].separation(coordExo).arcsecond
      ind = np.where(sep <= np.nanmin(sep))[0]
      if len(ind) > 0:
        outres = exo.loc[ind]
    else:
      sep2 = coordSCNASA[SCind].separation(coordNasa).arcsecond
      ind2 = np.where(sep2 <= np.nanmin(sep2))[0]
      if len(ind2) > 0:
        outres = nasa_eu.loc[ind2]
    return outres

  def get_mass_taken(row):
    if row['mass'] > 0:
      return row['mass']
    if row['mass_sini'] > 0:
      return row['mass_sini']
    return None

  def get_mass_taken_error(row):
    if row['mass'] > 0:
      return row['mass_error_max']
    if row['mass_sini'] > 0:
      return row['mass_sini_error_max']
    return None


  EXOcsv = "../exo.csv"
  NASAcsv = "../nasaexo.csv"
  SCcsv = "download/SWEETCAT_Dataframe.csv"
  dtype_SW = {'gaia_dr2':'int64','gaia_dr3':'int64'}

  exo = pd.read_csv(EXOcsv)
  exo = exo[(exo.detection_type == 'Radial Velocity') | (exo.detection_type == 'Primary Transit') | (exo.detection_type == 'Astrometry')]
  exo = exo.reset_index()
  nasa = pd.read_csv(NASAcsv)
  SC = pd.read_csv(SCcsv, dtype=dtype_SW)
  print(f"lines in exo: {len(exo)}")
  print(f"lines in NASA: {len(nasa)}")
  print(f"lines in SC: {len(SC)}")

  # Change columns' names to be equal to EU (more readable names)
  old_col=[val for key,val in dic.items() if val!='']
  new_col =[key for key,val in dic.items() if val!='']
  nasa_eu = nasa.rename(columns=dict(zip(old_col, new_col)))

  #Getting coordinates matching
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
  coordNasa   = coord.SkyCoord(ra=nasa_eu['ra'].values,
                             dec=nasa_eu['dec'].values,
                             unit=(u.deg, u.deg),
                             frame='icrs')

  #create planets dataframe
  if load:
    data = pd.read_csv("data_new.csv")
  else:
    column_names = ["name","swname", "mass", "mass_error_max","mass_sini", "mass_sini_error_max", "mass_taken", "mass_taken_error", "orbital_period", "orbital_period_error_max", "orbital_period_error_min", "swfeh", "swefeh", "swteff", "pl_source","swsource", "swMass"]
    data = pd.DataFrame([], columns =column_names)
    for i in range(len(SC)):
      planets = get_planets_for_SC_star(i)
      planets["mass_taken"] = planets.apply(lambda row: get_mass_taken(row), axis=1)
      planets["mass_taken_error"] = planets.apply(lambda row: get_mass_taken_error(row), axis=1)
      planets["swname"] = SC.Name[i]
      planets["swfeh"] = SC["[Fe/H]"][i]
      planets["swefeh"] = SC["e[Fe/H]"][i]
      planets["swteff"] = SC.Teff[i]
      planets["pl_source"] = "EU"
      planets["swsource"] = SC.SWFlag[i]
      planets["swMass"] = SC.Mass_t[i]
      dtmp = planets[column_names]
      data = data.append(dtmp, ignore_index=True)
    data.to_csv("data_new.csv", index=False, header=True)





### Main program:
def main():


  #sweet_database = "/home/sousasag/Nextcloud/WORK/spectra/sweet_cat_spec_2020/prog/database/SWEETCAT_Dataframe_fix.csv"
  sweet_database = "download/SWEETCAT_Dataframe.csv"
  SC_full_load = pd.read_csv(sweet_database, dtype=dtype_SW)
  #create_html_table(SC_full_load)
  #SC_stats(SC_full_load, plot_stats=True)
  #plot_sky_SC_chart(SC_full_load)

  #plot_histograms_parameters(SC_full_load)

  #see in the jupyter notebook (SWEET-CAT_Exoeu_link.ipynb) how to create this file
  #get_planet_propreties(load = True)
  plot_feh_mass_period_all("data_new.csv")
  return





if __name__ == "__main__":
    main()