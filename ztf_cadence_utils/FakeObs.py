import pandas as pd
from pandas import DataFrame
import numpy as np
from random import randint
from random import seed

def generateObsData(cad_g, cad_r, cad_i, skynoise_g, skynoise_r, skynoise_i, MJD_min, zp, healpixID_min, healpixID_max, step, seed):
    
    df = pd.DataFrame()
    df_r = pd.DataFrame()
    df_i = pd.DataFrame()
    df_tot = pd.DataFrame()

    MJD = np.arange(MJD_min, MJD_min+180, cad_g)
    healpixID = np.arange(healpixID_min, healpixID_max, step)
    
    for pix in healpixID:
        
        np.random.seed(seed)
    
        df['time'] = MJD 
        df['field'] = np.random.randint(250, 1866)
        df['skynoise'] = skynoise_g
        df['zp'] = zp
        df['rcid'] = 1.0
        df['band'] = 'ztfg'
        #df['ccd_x'] = 1.0
        df['programid'] = 1.0
        df['ra'] = np.random.randint(35, 72)
        df['dec'] = np.random.randint(-27, 79)
        df['i'] = -1.0
        df['j'] = -1.0
        df['ra_quad'] = df['ra']
        df['dec_quad'] = df['dec']
        #df['ccd_y'] = 1.0
        df['healpixID'] = pix
        df['healpixID']= df['healpixID'].map(str)

        MJD_r = np.arange(MJD_min+(15/(25*60)), MJD_min+(15/(25*60))+180, cad_r)

        df_r['time'] = MJD_r
        df_r['field'] = np.random.randint(250, 1866)
        df_r['skynoise'] = skynoise_r
        df_r['zp'] = zp
        df_r['rcid'] = 1.0
        df_r['band'] = 'ztfr'
        #df_r['ccd_x'] = 1.0
        df_r['programid'] = 1.0
        df_r['ra'] = np.random.randint(35, 72)
        df_r['dec'] = np.random.randint(-27, 79)
        df_r['i'] = -1.0
        df_r['j'] = -1.0
        df_r['ra_quad'] = df_r['ra']
        df_r['dec_quad'] = df_r['dec']
        #df_r['ccd_y'] = 1.0
        df_r['healpixID'] = pix
        df_r['healpixID']= df_r['healpixID'].map(str)

        MJD_i = np.arange(MJD_min+(30/(25*60)), MJD_min+(30/(25*60))+180, cad_i)

        df_i['time'] = MJD_i
        df_i['field'] = np.random.randint(250, 1866)
        df_i['skynoise'] = skynoise_i
        df_i['zp'] = zp
        df_i['rcid'] = 1.0
        df_i['band'] = 'ztfi'
        #df_i['ccd_x'] = 1.0
        df_i['programid'] = 1.0
        df_i['ra'] = np.random.randint(35, 72)
        df_i['dec'] = np.random.randint(-27, 79)
        df_i['i'] = -1.0
        df_i['j'] = -1.0
        df_i['ra_quad'] = df_i['ra']
        df_i['dec_quad'] = df_i['dec']
        #df_i['ccd_y'] = 1.0
        df_i['healpixID'] = pix
        df_i['healpixID']= df_i['healpixID'].map(str)

        df_tot = pd.concat([df_tot, df, df_r, df_i])
    
    return df_tot
    