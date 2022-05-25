import numpy as np
import astropy.units as u
import pandas as pd
import healpy as hp
from ztf_metrics.metricUtils import dustMap, seasons, addNight, coaddNight
from ztf_pipeutils.ztf_hdf5 import Read_LightCurve
from ztf_simfit_plot.z_bins import Z_bins


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

    def run(self, pixnum, data):
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
        dfm = data.groupby(['season']).apply(
            lambda x: pd.DataFrame({'time_min': [x['time'].min()], 'time_max': [x['time'].max()]})).reset_index()
        df = df.merge(dfm, left_on=['season'], right_on=['season'])
        # df['time_max'] = data.groupby(['season'])['time'].max()

        print(data.groupby(
            ['season']).apply(lambda x: x['time'].min()))
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


class RedMagMetric:
    def __init__(self, gap=60, nside=128, coadd_night=1,
                 meta_dir='/Users/manon/ztf_pixelized/meta_pix', plot=False):

        self.gapvalue = gap
        self.nside = nside
        self.coadd_night = coadd_night
        self.dustmap = dustMap(nside)
        self.meta_dir = meta_dir
        self.plot = plot

    def run(self, pixnum, data):

        self.pixnum = pixnum

        mk_sup = data['c_err'] > 0.036
        mk_inf = data['c_err'] < 0.044
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
        df_new = pd.DataFrame()

        dataFileName = 'Meta_fit_{}_{}.hdf5'.format(pixnum, s[0])

        m, delta_m = self.get_mag(data['z'])
        mag = np.mean(m)
        delta_mag = np.mean(delta_m)

        df_ = self.get_df(dataFileName)
        df_['healpixID'] = pixnum
        df_['mag_lim'] = mag
        df_['mag_lim_max'] = mag + delta_mag
        df_['mag_lim_min'] = mag - delta_mag

        df_new = df.merge(df_, on='healpixID')

        return df_new

    def get_dl(self, z):

        from astropy.cosmology import FlatLambdaCDM

        cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc,
                              Tcmb0=2.725 * u.K, Om0=0.3)
        dl_Mpc = cosmo.luminosity_distance([z])
        dl_pc = dl_Mpc.to(u.pc)

        h = cosmo.h
        delta_dl_Mpc = (cosmo.luminosity_distance(
            [z]+h) - cosmo.luminosity_distance([z]-h))/2*h
        delta_dl_pc = delta_dl_Mpc.to(u.pc)

        return dl_pc, delta_dl_pc

    def get_mag(self, z):

        dl_pc, delta_dl_pc = self.get_dl(z)

        m = -2.5*np.log10(dl_pc.to_value()**2/10**2) + 19.6
        delta_m = delta_dl_pc / dl_pc

        return m, delta_m

    def get_df(self, file_name_):

        z_comp, z_comp_sup, z_comp_inf = self.get_z(self.meta_dir, file_name_)

        d = {'z_comp': [z_comp], 'z_comp_max': [
            z_comp_sup], 'z_comp_min': [z_comp_inf]}
        return pd.DataFrame(d)

    def get_z(self, diretory, filename):

        z_comp, z_comp_sup, z_comp_inf = Z_bins(metaFitInput=filename,
                                                inputDir=diretory).get_z()
        if self.plot:
            Z_bins(metaFitInput=filename, inputDir=diretory).plot_z_comp()

        return z_comp, z_comp_sup, z_comp_inf
