import pandas as pd
from pandas import DataFrame
import numpy as np

def generateObsData(cad_g, cad_r, cad_i, skynoise_g, skynoise_r, skynoise_i, MJD_min, zp, healpixID_min, healpixID_max, step):
    
    df = pd.DataFrame()
    MJD = np.arange(MJD_min, MJD_min+180, cad_g)
    print(MJD, len(MJD))
    healpixID = np.arange(healpixID_min, len(MJD), step)
    print(healpixID, len(healpixID))
    
    df['time'] = MJD 
    df['field'] = 1.0
    df['skynoise'] = skynoise_g
    df['zp'] = zp
    df['rcid'] = 1.0
    df['band'] = 'ztfg'
    df['ccd_x'] = 1.0
    df['programid'] = 1.0
    df['ra'] = 71
    df['dec'] = -16
    df['i'] = -1.0
    df['j'] = -1.0
    df['ra_quad'] = 69
    df['dec_quad'] = -16
    df['ccd_y'] = 1.0
    df['healpixID'] = healpixID
    
    df_r = pd.DataFrame()
    MJD_r = np.arange(MJD_min+(15/(25*60)), MJD_min+(15/(25*60))+180, cad_r)
    healpixID_r = np.arange(healpixID_min, len(MJD_r), step)    
    
    df_r['time'] = MJD_r
    df_r['field'] = 1.0
    df_r['skynoise'] = skynoise_r
    df_r['zp'] = zp
    df_r['rcid'] = 1.0
    df_r['band'] = 'ztfr'
    df_r['ccd_x'] = 1.0
    df_r['programid'] = 1.0
    df_r['ra'] = 71
    df_r['dec'] = -16
    df_r['i'] = -1.0
    df_r['j'] = -1.0
    df_r['ra_quad'] = 69
    df_r['dec_quad'] = -16
    df_r['ccd_y'] = 1.0
    df_r['healpixID'] = healpixID_r
    
    df_i = pd.DataFrame()
    MJD_i = np.arange(MJD_min+(30/(25*60)), MJD_min+(30/(25*60))+180, cad_i)
    healpixID_i = np.arange(healpixID_min, len(MJD_i), step)
    
    df_i['time'] = MJD_i
    df_i['field'] = 1.0
    df_i['skynoise'] = skynoise_i
    df_i['zp'] = zp
    df_i['rcid'] = 1.0
    df_i['band'] = 'ztfi'
    df_i['ccd_x'] = 1.0
    df_i['programid'] = 1.0
    df_i['ra'] = 71
    df_i['dec'] = -16
    df_i['i'] = -1.0
    df_i['j'] = -1.0
    df_i['ra_quad'] = 69
    df_i['dec_quad'] = -16
    df_i['ccd_y'] = 1.0
    df_i['healpixID'] = healpixID_i
    
    df_tot = pd.concat([df, df_r, df_i])
    
    return df_tot
    