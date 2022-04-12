import pandas as pd
from optparse import OptionParser
import numpy as np

class CadenceMetric:
    def __init__(self, gap=60):
        self.gap_ = gap
        
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
        index = list((df_time[df_time>self.gap_].index)-1)
        index.insert(0, df.index.min())
        index.insert(len(index), df.index.max())
        df['season'] = 0
        
        for i in range(0, len(index)-1):
            df.loc[index[i] : index[i+1], 'season'] = i+1    
        return df
    
    def run(self, data, numpix):
        """
        Running method.
        
        Parameters
        --------------
        df: pandas df
            data to process with observations.
        numpix : int
            numero of the pixel related to the observation.
        
        Returns
        --------------
        A new data frame with the desired parameters. It has as many lines as the observation season for the selected pixel.
        """
        data = self.season_l(data)
        s = data['season'].unique()
        df = pd.DataFrame(s, columns=['season'])
        df['healpixID'] = numpix
        print('DF=', df)
               
        df_b = data.groupby(['season']).apply(lambda x : self.calc_metric(group=x))
        df = df.merge(df_b, left_on=['season'], right_on=['season'])
        for b in ['ztfg', 'ztfr', 'ztfi']:
            df_b = data.groupby(['season']).apply(lambda x : self.calc_metric(group=x,bands=[b],colnames=
                               ('cad_{}'.format(b),'nb_obs_{}'.format(b), 'gap_{}'.format(b),
                               'season_lenght_{}'.format(b), 'skynoise_{}'.format(b)), 
                               to_calc=['cadence', 'nb_obs', 'gap', 'season_lenght', 'skynoise']))
            df = df.merge(df_b, left_on=['season'], right_on=['season'])
            
        return df
            
    def calc_metric(self, group, bands=['ztfg', 'ztfr', 'ztfi'], colnames=('cad_all', 'nb_obs_all', 
                             'gap_all', 'season_lenght_all', 'skynoise_all'), 
                                                to_calc = ['cadence', 'nb_obs', 'gap', 'season_lenght', 'skynoise']):
        
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
            self.grp = grp
            
            for calc in to_calc:
                res[corresp[calc]]=[eval('self.{}()'.format(calc))]
        return pd.DataFrame.from_dict(res)
    
    def cadence(self):
        """
        Calculation of the candence.
       
        Returns
        --------------
        Value of the candence.
        """
        diff = self.grp['time'].diff()
        cad = diff.median()
        return cad
    
    def nb_obs(self):
        """
        Calculates the number of observations.
        
        Returns
        --------------
        Value of the number of observations.
        """
        nb = len(self.grp)
        return nb
    
    def gap(self):
        """
        Calculation of the candence.
        
        Returns
        --------------
        Value of the season lenght gap.
        """
        diff = self.grp['time'].diff()
        g = diff.max()
        return g
    
    def season_lenght(self):
        """
        Calculates the season lenght.
        
        Returns
        --------------
        Value of the season lenght.
        """
        
        sl = self.grp['time'].max() - self.grp['time'].min()
        return sl
    
    def skynoise(self):
        """
        Calculates the skynoise.
        
        Returns
        --------------
        Value of the skynoise mean.
        """
        
        sk = self.grp['skynoise'].mean()
        return sk


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
    for b in bo_ok[:10]:
        df_ = df[df['healpixID'].str.contains(b)]
        df_new = df_.copy()
        pd_ = cl.run(df_new, b)

        df1 = pd.concat([df1, pd_])
    return df1
        
        
        
        
        
        
        
        
        
        
        



