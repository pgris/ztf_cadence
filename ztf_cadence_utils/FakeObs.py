import pandas as pd
from pandas import DataFrame
import numpy as np
from random import randint
from random import seed


def generateObsData(cad_g, cad_r, cad_i, N_g, N_r, N_i, skynoise_g, skynoise_r, skynoise_i, MJD_min, zp_g, zp_r, zp_i, healpixID_min,
                    healpixID_max, step, seed):

    df_tot = pd.DataFrame()

    MJD = [np.arange(MJD_min, MJD_min+180, cad_g), np.arange(MJD_min+(15/(25*60)), MJD_min+(15/(25*60))+180, cad_r),
           np.arange(MJD_min+(30/(25*60)), MJD_min+(30/(25*60))+180, cad_i)]
    healpixID = 70905
    skynoise = [skynoise_g, skynoise_r, skynoise_i]
    nvisits = [N_g, N_r, N_i]
    zp = [zp_g, zp_r, zp_i]
    band = ['ztfg', 'ztfr', 'ztfi']

    for i, band in enumerate(band):
        df = pd.DataFrame()
        print(band)

        df['time'] = MJD[i]
        df['field'] = 401
        df['skynoise'] = skynoise[i]+1.25*np.log10(nvisits[i])
        df['zp'] = zp[i]
        df['rcid'] = 45
        df['band'] = band
        df['programid'] = 1.0
        df['ra'] = 33.398438
        df['dec'] = -1.193748
        df['i'] = 0
        df['j'] = -1.0
        df['ra_quad'] = 33.706826
        df['dec_quad'] = -1.377433
        df['healpixID'] = healpixID
        df['healpixID'] = df['healpixID'].map(str)

        df_tot = pd.concat([df_tot, df])

    return df_tot
