import pandas as pd
import healpy as hp
import numpy as np
from ztf_metrics.metricUtils import dustMap, seasons, addNight, coaddNight
from ztf_pipeutils.ztf_pipeutils.ztf_hdf5 import Read_LightCurve
from ztf_simfit.ztf_simfit_plot.z_bins import Z_bins


class CadenceMetric:
    def __init__(self, gap=60, nside=128, coadd_night=1):
        """
        Class that allows to compute different observables like cadence and season length from a pandas data frame.

        Parameters
        --------------
        gap: int, opt
            gap between two seasons (default : 60 days)
        nside: int,opt
          healpix nside parameter (default: 128)
        coadd_night: int, opt
          to perform coaddition per band and per night

        """
        self.gapvalue = gap
        self.nside = nside
        self.coadd_night = coadd_night
        self.dustmap = dustMap(nside)

    def run(self, pixnum, data, input_dir_data):
        """
        Running method.

        Parameters
        --------------
        pixnum : int
            pixel ID related to the observation.
         df: pandas df
            data to process with observations.

        Returns
        --------------
        A new data frame with the desired parameters. It has as many lines as the observation season for the selected pixel.
        """

        # modify healpixID
        data['healpixID'] = pixnum
        # define nights
        data = addNight(data)
        # coadd here
        if self.coadd_night:
            data = coaddNight(data, cols=['night', 'band', 'healpixID'])

        # get E(B-V)
        idx = self.dustmap['healpixID'] == pixnum
        seldust = self.dustmap[idx]

        # get seasons
        data = seasons(data, self.gapvalue, mjdCol='time')
        s = data['season'].unique()
        df = pd.DataFrame(s, columns=['season'])
        df['season'] = df['season'].astype(int)
        df['healpixID'] = pixnum
        df['ebvofMW'] = seldust['ebvofMW'].to_list()[0]
        df['healpixRA'] = seldust['RA'].to_list()[0]
        df['healpixDec'] = seldust['Dec'].to_list()[0]

        # need to coadd by night here
        data_coadded = coaddNight(data.drop(columns=['band']), cols=[
                                  'night', 'healpixID'])
        data_coadded['band'] = 'ztfall'
        df_b = data_coadded.groupby(['season']).apply(
            lambda x: self.calc_metric(group=x, bands=['ztfall'])).reset_index()
        df = df.merge(df_b, left_on=['season'], right_on=['season'])
        for b in ['ztfg', 'ztfr', 'ztfi']:
            df_b = data.groupby(['season']).apply(lambda x: self.calc_metric(group=x, bands=[b], colnames=('cad_{}'.format(b), 'nb_obs_{}'.format(b), 'gap_{}'.format(b),
                                                                                                           'season_length_{}'.format(b), 'skynoise_{}'.format(b)),
                                                                             to_calc=['cadence', 'nb_obs', 'gap', 'season_length', 'skynoise']))
            df = df.merge(df_b, left_on=['season'], right_on=['season'])

        return df

    def calc_metric(self, group, bands=['ztfg', 'ztfr', 'ztfi'],
                    colnames=('cad_all', 'nb_obs_all',
                              'gap_all', 'season_length_all', 'skynoise_all'),
                    to_calc=['cadence', 'nb_obs', 'gap', 'season_length', 'skynoise']):
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
                             'gap_all', 'season_length_all'))
        to_calc: list, opt
            list of the different parameters you want to calculate (default : ['cadence', 'nb_obs', 'gap', 'season_length']).
            the parameters you can put in this list : 'cadence', 'nb_obs', 'gap', 'season_length', 'skynoise'.
            the length of to_talc have to match with the length of colnames.

        Returns
        --------------
        A new data frame with the desired calculates parameters.
        """

        idx = group['band'].isin(bands)
        grp = group[idx]

        corresp = dict(zip(to_calc, colnames))
        res = {}

        if len(grp) <= 1:
            for calc in to_calc:
                res[corresp[calc]] = [0]
        else:
            grp.sort_values('time', ascending=True)
            self.grp = grp

            for calc in to_calc:
                res[corresp[calc]] = [eval('self.{}()'.format(calc))]

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
        Value of the season length gap.
        """
        diff = self.grp['time'].diff()
        g = diff.max()
        return g

    def season_length(self):
        """
        Calculates the season length.

        Returns
        --------------
        Value of the season length.
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

    
class RedshiftCompMetric:
    def __init__(self, gap=60, nside=128, coadd_night=1):
        
        self.gapvalue = gap
        self.nside = nside
        self.coadd_night = coadd_night
        self.dustmap = dustMap(nside)

    def run(self, pixnum, data, input_dir_data):
        
        self.input_dir_data = input_dir_data
        self.pixnum = pixnum
        
        mk_sup = data['c_err']>0.036
        mk_inf = data['c_err']<0.044
        mask_c_err = mk_sup & mk_inf
        
        # get E(B-V)
        idx = self.dustmap['healpixID'] == pixnum   
        seldust = self.dustmap[idx]
        
        # get seasons
        s = pd.unique(data['season'])
        df = pd.DataFrame(s, columns=['season'])
        df['season'] = df['season'].astype(int)
        df['healpixID'] = pixnum
        
        data = data[mask_c_err]
        
        dataFileName = data.meta['file_name']
        lc = self.get_lc(dataFileName, data)
        df_new = pd.DataFrame()
        
        for b in pd.unique(lc['band']):
            df_b = self.get_df(lc, data.meta['directory'], data.meta['file_name_meta'], b)
            df_b['healpixID'] = pixnum
            
            df_c = df.merge(df_b, on ='healpixID')
            
            df_new = pd.concat([df_new, df_c], ignore_index=True)
       
        return df_new
            
        
    def get_df(self, lc, dir_, file_name_, b):
        
        dico={'variables': {'sel': 'sel','c_err': 'c_err','zmin': 'z','zmax': 'z','chi2': 'chi2'},
              'operators': {'op sel': 'operator.eq','op c_err': 'operator.ne','op zmin': 'operator.ge',
                            'op zmax': 'operator.lt','op_chi2': 'operator.lt'},
              'limites': {'lim sel': 1,'lim c_err': -1,'lim zmin': 0.01,'lim zmax': 0.1,'lim_chi2': 2}}
        
        mag_min = self.calc_mag(lc, b)
        z_comp, z_comp_sup, z_comp_inf = self.get_z(dir_, file_name_, dico)
        
        return pd.DataFrame({'z_comp' : z_comp, 'z_comp_sup' : z_comp_sup, 'z_comp_inf' : z_comp_inf,
                            'mag_min' : mag_min, 'band' : b})
            

    def get_lc(self, dataFileName, data):
        cl2 = Read_LightCurve(file_name=dataFileName, inputDir=self.input_dir_data)
        for path in data['path']:
            lc = cl2.get_table(path=path)
            lc = self.complete_lc(lc, 5)
            
            return lc
    
    def complete_lc(self, lc, snr):

        lc['SNR'] = lc['flux'] / lc['fluxerr']
        lc['phase'] = (lc['time'] - lc.meta['t0']) / (1-lc.meta['z'])
        idx = lc['SNR'] >= snr

        return lc[idx]   
        
    def calc_mag(self, lc, band):
        
        mask_b = lc['band'] == band
        lc = lc[mask_b]
        flux_max = np.max(lc['flux'])
        mag_min = -2.5*np.log10(flux_max) + lc[lc['flux']==flux_max]['zp']

        return mag_min   
    
    def get_z(self, diretory, filename, dic):
        
        z_comp, z_comp_sup, z_comp_inf = Z_bins(metaFitInput = filename, 
                                                inputDir = diretory, dico = dic).plot_interpolation1d(add_column=True)
        
        return z_comp, z_comp_sup, z_comp_inf
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        