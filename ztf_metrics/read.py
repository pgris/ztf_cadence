import pandas as pd
from optparse import OptionParser
import matplotlib.pyplot as plt
import ligo.skymap.plot
import numpy as np
from matplotlib.projections import get_projection_names
from astropy.coordinates import SkyCoord
from astropy import units as u
import healpy as hp
from itertools import chain

parser = OptionParser()

parser.add_option('--fileName', type=str, default='data_36.0_72.0.hdf5',
                  help='meta data file name [%default]')
parser.add_option('--input_dir', type=str, default='ztf_pixelized',
                  help='folder directory name [%default]')


opts, args = parser.parse_args()

fileName = opts.fileName
input_dir = opts.input_dir

class CadenceMetric:
    def __init__(self, SL=60):
        self.SL = SL
    
    def season_l(self, df):
        df = df.sort_values('time', ascending=True)
        df = df.reset_index()
        df_time = df['time'].diff()
        index = list((df_time[df_time>self.SL].index)-1)
        index.insert(0, df.index.min())
        index.insert(len(index), df.index.max())
        #print('df_time max=', np.max(df_time), 'index=', index)
        df['season'] = 0
        
        for i in range(0, len(index)-1):
            df.loc[index[i] : index[i+1], 'season'] = i+1    
        return df
        
    def run(self, data):
        data = self.season_l(data)
        s = data['season'].unique()
        df = pd.DataFrame(s, columns=['season'])
               
        df_b = data.groupby(['season']).apply(lambda x : self.cadence(group=x))
        df = df.merge(df_b, left_on=['season'], right_on=['season'])
        for b in ['ztfg', 'ztfr', 'ztfi']:
            df_b = data.groupby(['season']).apply(lambda x : self.cadence(group=x,bands=[b],colnames=
                               ('cad_{}'.format(b),'nb_obs_{}'.format(b), 'gap_{}'.format(b), 'nb_obs_{}'.format(b),
                               'season_lenght_{}'.format(b))))
            df = df.merge(df_b, left_on=['season'], right_on=['season'])
            
        return df
            
    def cadence(self, group, bands=['ztfg', 'ztfr', 'ztfi'], colnames=('cad_all', 'nb_obs_all', 'gap_all', 'nb_obs_all', 
                                                                       'season_lenght_all')):
        idx = group['band'].isin(bands)
        grp = group[idx]
        if len(grp)<=1 :
            cd = -1
            SL = -1
            obs = len(grp)
            gap = -1
        else :
            grp.sort_values('time', ascending=True)
            diff = grp['time'].diff()
        
            cd = diff.median()
            SL = grp['time'].max() - grp['time'].min()
            obs = len(grp)
            gap = np.max(diff)
        
        return pd.DataFrame({colnames[0] : [cd], colnames[1]: [obs], colnames[2]: [gap], colnames[3]: [obs], colnames[4]: [SL]})
    
    
df = pd.read_hdf('{}/{}'.format(input_dir, fileName))
bo = ','.join(df['healpixID'].to_list())
bo_= bo.split(",")
bo_ok = [i for i in bo_ if i!='None']
bo_ok = sorted(list(set(bo_ok)))

cl = CadenceMetric()   

df1 = pd.DataFrame(columns=['healpixID','season','cad_all','nb_obs_all','gap_all','cad_ztfg','nb_obs_ztfg','gap_ztfg',
                           'cad_ztfr','nb_obs_ztfr','gap_ztfr','cad_ztfi','nb_obs_ztfi','gap_ztfi'])

for i in range(200,210):
    df_ = df[df['healpixID'].str.contains(bo_ok[i])]
    df_new = df_.copy()
    pd_ = cl.run(df_new)
    pd_.insert(loc = 0, 
          column = 'healpixID', 
          value = bo_ok[i]) 
    
    df1 = pd.concat([df1, pd_])

print(df1)

        
        
        
        
        
        
        
        
        
        
        
        
        


