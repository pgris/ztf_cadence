import pandas as pd
from optparse import OptionParser
import numpy as np

class CadenceMetric:
    def __init__(self, GAP=60):
        self.GAP = GAP
        
        """
        Class that allows to compute different observables like cadence and season length from a pandas data frame.
        
        Parameters
        --------------
        GAP: int, opt
            gap for the season length in days (default : 60)
  
        """
        
    def season_l(self, df):
        """
        Method that determines the different seasons of an observation related to a pixel.
        
        Parameters
        --------------
        df: pandas df
            data to process with observations.
        
        Returns
        --------------
        Starting pandas with a 'season' column.
        """
        
        df = df.sort_values('time', ascending=True)
        df = df.reset_index()
        df_time = df['time'].diff()
        index = list((df_time[df_time>self.GAP].index)-1)
        index.insert(0, df.index.min())
        index.insert(len(index), df.index.max())
        df['season'] = 0
        
        for i in range(0, len(index)-1):
            df.loc[index[i] : index[i+1], 'season'] = i+1    
        return df
    
    def run(self, data):
        """
        Running method.
        
        Parameters
        --------------
        df: pandas df
            data to process with observations.
        
        Returns
        --------------
        A new data frame with the desired parameters. It has as many lines as the observation season for the selected pixel.
        """
        data = self.season_l(data)
        s = data['season'].unique()
        df = pd.DataFrame(s, columns=['season'])
               
        df_b = data.groupby(['season']).apply(lambda x : self.calc_metric(group=x))
        df = df.merge(df_b, left_on=['season'], right_on=['season'])
        for b in ['ztfg', 'ztfr', 'ztfi']:
            df_b = data.groupby(['season']).apply(lambda x : self.calc_metric(group=x,bands=[b],colnames=
                               ('cad_{}'.format(b),'nb_obs_{}'.format(b), 'gap_{}'.format(b),
                               'season_lenght_{}'.format(b)), to_calc=['cadence', 'nb_obs', 'gap', 'season_lenght']))
            df = df.merge(df_b, left_on=['season'], right_on=['season'])
            
        return df
            
    def calc_metric(self, group, bands=['ztfg', 'ztfr', 'ztfi'], colnames=('cad_all', 'nb_obs_all', 
                             'gap_all', 'season_lenght_all'), to_calc = ['cadence', 'nb_obs', 'gap', 'season_lenght']):
        
        """
        Method that calculates.
        
        Parameters
        --------------
        group: pandas df
            data frame group by season.
        bands: list, opt
            list of the different bands in your data frame (default : ['ztfg', 'ztfr', 'ztfi'])
        colnames: list, opt
            list of the different colnames for your futur data frame (default : ('cad_all', 'nb_obs_all', 
                             'gap_all', 'season_lenght_all'))
        to_calc: list, opt
            list of the different parameters you want to calculate (default : ['cadence', 'nb_obs', 'gap', 'season_lenght']).
            the parameters you can put in this list : 'cadence', 'nb_obs', 'gap', 'season_lenght', 'skynoise'.
            the lenght of to_talc have to match with the lenght of colnames.
        
        Returns
        --------------
        A new data frame with the desired calculates parameters.
        """
        
        idx = group['band'].isin(bands)
        grp = group[idx]
        
        corresp = dict(zip(to_calc, colnames)) 
        res = {}
        
        if len(grp)<=1 :
            for calc in to_calc:
                res[corresp[calc]]=[0]
        else :
            grp.sort_values('time', ascending=True)
            diff = grp['time'].diff()
            
            obj_for_calc = {}
            obj_for_calc['cadence'] = diff
            obj_for_calc['nb_obs'] = grp
            obj_for_calc['gap'] = diff
            obj_for_calc['season_lenght'] = grp
            
            for calc in to_calc:
                res[corresp[calc]]=[eval('self.{}(obj_for_calc[\'{}\'])'.format(calc, calc))]
        return pd.DataFrame.from_dict(res)
    
    def cadence(self, diff):
        """
        Calculation of the candence.
        
        Parameters
        --------------
        diff: pandas 
            data corresponding to a difference in time.
        Returns
        --------------
        Value of the candence.
        """
        cad = diff.median()
        return cad
    
    def nb_obs(self, grp):
        """
        Calculates the number of observations.
        
        Parameters
        --------------
        grp: pandas 
            data corresponding to the data frame group by season and bands.
        Returns
        --------------
        Value of the number of observations.
        """
        nb = len(grp)
        return nb
    
    def gap(self, diff):
        """
        Calculation of the candence.
        
        Parameters
        --------------
        diff: pandas 
            data corresponding to a difference in time.
        Returns
        --------------
        Value of the season lenght gap.
        """
        g = diff.max()
        return g
    
    def season_lenght(self, grp):
        """
        Calculates the season lenght
        
        Parameters
        --------------
        grp: pandas 
            data corresponding to the data frame group by season and bands.
        Returns
        --------------
        Value of the season lenght.
        """
        
        sl = grp['time'].max() - grp['time'].min()
        return sl


def read(input_dir, fileName):
    """
    Method which read the pandas data frame of observations for all the pixels.

    Parameters
    --------------
    input_dir: str, opt 
        input directory.
    fileName: str, opt 
        file name.
    Returns
    --------------
    New data frame for the differents pixels with the calculation (value) of differents metrics.
    """
    
    df = pd.read_hdf('{}/{}'.format(input_dir, fileName))
    bo = ','.join(df['healpixID'].to_list())
    bo_= bo.split(",")
    if 'None' in bo_:
        bo_.remove('None')
        bo_ = set(bo_)
    bo_ok = list(bo_)

    cl = CadenceMetric()   

    df1 = pd.DataFrame()
    for b in bo_ok[10:20]:
        df_ = df[df['healpixID'].str.contains(b)]
        df_new = df_.copy()
        pd_ = cl.run(df_new)
        pd_.insert(loc = 0, column = 'healpixID', value = b) 

        df1 = pd.concat([df1, pd_])
    return df1
        
        
        
        
        
        
        
        
        
        
        



