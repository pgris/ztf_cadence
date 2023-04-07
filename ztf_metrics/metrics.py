import numpy as np
import astropy.units as u
import pandas as pd
import healpy as hp
from ztf_metrics.metricUtils import dustMap, seasons, addNight, coaddNight
from ztf_pipeutils.ztf_hdf5 import Read_LightCurve
from ztf_simfit_plot.z_bins import Z_bins
from ztf_metrics.plotUtils import binnedData
from scipy import interpolate


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
        dfm = data.groupby('season').apply(
            lambda x: pd.DataFrame({'time_min': [x['time'].min()],
                                    'time_max': [x['time'].max()]})).reset_index()

        df = df.merge(dfm, left_on=['season'], right_on=['season'])
        # df['time_max'] = data.groupby(['season'])['time'].max()

        # need to coadd by night here
        data_coadded = coaddNight(data.drop(columns=['band']), cols=[
            'night', 'healpixID'])
        data_coadded['band'] = 'ztfall'
        data_coadded['season'] = data_coadded['season'].astype(int)
        df_b = data_coadded.groupby(['season']).apply(
            lambda x: self.calc_metric(group=x, bands=['ztfall'])).reset_index()
        df = df.merge(df_b, left_on=['season'], right_on=['season'])

        to_calc = ['cadence', 'nb_obs', 'gap', 'season_length', 'skynoise']
        for b in ['ztfg', 'ztfr', 'ztfi']:
            colnames = ('cad_{}'.format(b), 'nb_obs_{}'.format(b),
                        'gap_{}'.format(b),
                        'season_length_{}'.format(b), 'skynoise_{}'.format(b))
            df_b = data.groupby(['season']).apply(lambda x:
                                                  self.calc_metric(group=x,
                                                                   bands=[b],
                                                                   colnames=colnames,
                                                                   to_calc=to_calc))
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
    def __init__(self, gap=60, nside=128, coadd_night=1):
        """
        class to estimate the redshift and magnitude limits

        Parameters
        --------------

        """

        self.gapvalue = gap
        self.nside = nside
        self.coadd_night = coadd_night
        self.dustmap = dustMap(nside)
        self.sigmaC = 0.04
        self.absPeakMag = -19.6

        # cosmology
        from astropy.cosmology import FlatLambdaCDM

        self.cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc,
                                   Tcmb0=2.725 * u.K, Om0=0.3)

    def run(self, pixnum, data, plot=False):
        """
        Main metric method

        Parameters
        ---------------
        pixnum: int
           pixel number
        data: pandas df
           data to process
        plot: bool, opt
          to plot sigmaC vs z (zlim estimation) (default: False)

        Returns
        ----------
        pandas df with the following columns
        zlim, zlim_m,zlim_p,zlim_b
        mag, mag_m, mag_p
        healpixID
        season
        """

        self.pixnum = pixnum

        idx = data['healpixID'] == 'p{}p'.format(pixnum)
        data = data[idx]

        seasons = data['season'].unique()
        ebvofMW = np.mean(data['mwebv'])

        resdf = pd.DataFrame()
        for seas in seasons:
            idb = data['season'] == seas
            data_season = data[idb]
            ddf = self.process_season(data_season, plot=plot)
            ddf['healpixID'] = pixnum
            ddf['season'] = seas
            ddf['ebvofMW'] = ebvofMW
            resdf = pd.concat((resdf, ddf))

        return resdf

    def process_season(self, data_season, plot=False):
        """
        method to estimate metrics per pixel and season

        Parameters
        --------------
        data_season: pandas df
           data for the season
        plot: bool, opt
          to plot sigmaC vs z (zlim estimation) (default: False)

        Returns
        ----------
        pandas df with the following columns
        zlim, zlim_m,zlim_p,zlim_b
        mag, mag_m, mag_p

        """

        zlims = self.zlim(data_season, plot=plot)
        mag = self.get_mag(zlims)
        zlims.update(mag)

        zzlim = {}
        for key, vals in zlims.items():
            zzlim[key] = [vals.item()]

        return pd.DataFrame.from_dict(zzlim)

    def zlim(self, data, plot=False):
        """
        Method to estimate redshift limit values

        Parameters
        --------------
        data: pandas df
           data to process
        plot: bool, opt
          to plot sigmaC vs z (default: False)

        Returns
        ----------
        dict with the following keys: zlim, zlim_p, zlim_m, zlim_b

        """
        idx = data['fitstatus'] == 'fitok'
        seldata = data[idx]

        if len(seldata) < 5:
            return dict(zip(['zlim', 'zlim_p', 'zlim_m', 'zlim_b'], [np.asarray(-1.0)]*4))

        zmin, zmax = seldata['z'].min(), seldata['z'].max()
        bins = np.arange(zmin, zmax, 0.01)
        bin_values = binnedData(seldata, bins, 'z', 'c_err')

        if len(bin_values) < 2:
            return dict(zip(['zlim', 'zlim_p', 'zlim_m', 'zlim_b'], [np.asarray(-1.0)]*4))

        bin_values['c_err_p'] = bin_values['c_err']+bin_values['c_err_std']
        bin_values['c_err_m'] = bin_values['c_err']-bin_values['c_err_std']
        # bin_values['z'] = bin_centers

        zlim = {}
        for tt in ['', '_p', '_m']:
            zlim['zlim{}'.format(tt)] = self.zlim_interp(
                bin_values['c_err{}'.format(tt)], bin_values['z'])

        zlim['zlim_b'] = np.array(0)
        idx = bin_values['c_err'] >= self.sigmaC
        if len(bin_values[idx]) == 0:
            zlim['zlim_b'] = np.array(1)

        if plot:
            self.plot_sigmaC_z(bin_values, zlim)

        return zlim

    def plot_sigmaC_z(self, bin_values, zlim):
        """
        Method to plot sigmaC vs z
        Parameters
        --------------
        bin_values: pandas df
          data to plot
        zlim: dict
          zlim values to draw on the plot

        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        # ax.errorbar(
        #    bin_values['z'], bin_values['c_err'], yerr=bin_values['c_err_std'])
        ax.plot(bin_values['z'], bin_values['c_err'], color='r')
        zmin = bin_values['z'].min()
        zmax = bin_values['z'].max()
        ax.plot([zmin, zmax], [self.sigmaC]*2, ls='dotted', color='k')
        ax.fill_between(bin_values['z'], bin_values['c_err_p'],
                        bin_values['c_err_m'], color='yellow')
        plt.show()

    def zlim_interp(self, x, y):
        """
        Method to interpolate data

        Parameters
        ---------------
        x: array
          x-axis values
        y: array
          y-axis values

        Returns
        ----------
        y-axis value corresponding to self.sigmaC

        """
        interp = interpolate.interp1d(x, y, bounds_error=False, fill_value=-1.)

        return interp(self.sigmaC)

    def get_mag(self, zlim):
        """
        Method to estimate apparent magnitudes

        Parameters
        ---------------
        zlim: dict
          dict of redshift limits

        Returns
        ----------
        dict with the folling keys: mag, mag_p,mag_m

        """
        dl = {}
        for tt in ['', '_p', '_m']:
            z = zlim['zlim{}'.format(tt)]
            if z > 0.:
                dl_Mpc = self.cosmo.luminosity_distance([z])
                dl_pc = dl_Mpc.to(u.pc)
                m = 5.*np.log10(dl_pc.to_value()/10) + self.absPeakMag
            else:
                m = -999.
            dl['mag{}'.format(tt)] = np.array(np.mean(m))

        return dl
