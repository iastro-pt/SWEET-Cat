#!/usr/bin/python
## My first python code

##imports:

import pandas as pd
import numpy as np

## My functions:



### Main program:
def main():
  print("Hello")
  dtype_SW = dtype={'gaia_dr2':'int64','gaia_dr3':'int64'}
  SC = pd.read_csv("download/SWEETCAT_Dataframe.csv", dtype=dtype_SW)

  columns_html = ['Name','RA','DEC','gaia_dr3','Gmag','Plx','Distance','Teff','eTeff','Logg','eLogg','[Fe/H]','e[Fe/H]','Vt','eVt','Logg_gaia','eLogg_gaia','Mass_t','eMass_t','Radius_t','eRadius_t','SWFlag','Reference', 'simbad_link', 'reference_link']
  columns_format_name = ['Gmag'   ,'Plx'    ,'Distance','Teff'   ,'eTeff'  ,'Logg'   ,'eLogg'  ,'[Fe/H]' ,'e[Fe/H]','Vt'     ,'eVt'    ,'Logg_gaia','eLogg_gaia','Mass_t' ,'eMass_t','Radius_t','eRadius_t']
  columns_format      = ['{:5.2f}','{:5.2f}','{:6.1f}' ,'{:5.0f}','{:5.0f}','{:4.2f}','{:4.2f}','{:4.2f}','{:4.2f}','{:4.2f}','{:4.2f}','{:5.3f}'  ,'{:5.3f}'   ,'{:5.3f}','{:5.3f}','{:5.3f}' ,'{:5.3f}'  ]



  SC2 = SC.copy()
  SC2 = SC2.rename(columns={'Link': 'reference_link'})
  slinks = ["http://simbad.u-strasbg.fr/simbad/sim-id?Ident=GAIA+DR2+"+str(gdr2)+"&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id" for gdr2 in SC2['gaia_dr2']]
  SC2['simbad_link'] = slinks
  for i, c in enumerate(columns_format_name):
    SC2[c] = SC[c].map(columns_format[i].format).replace("   nan","").replace("  nan","").replace(" nan", "").replace("nan","") 

  SC2.to_csv("download/SWEETCAT_Dataframe_html.csv", index=False, columns = columns_html)

if __name__ == "__main__":
    main()

