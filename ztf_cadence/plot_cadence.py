import matplotlib.pyplot as plt
import healpy as hp
import matplotlib as mpl
import numpy as np
from ztf_pipeutils.ztf_util import checkDir, multiproc
import pandas as pd


class PlotOS:
    """
    class to plot OS pointings vs time

    Parameters
    --------------
    nside: int, opt
      nside for healpix (default: 128)
    outDir: str,opt
      output directory (default: plot_cadence)

    """

    def __init__(self, nside=128, outDir='plot_cadence'):

        # get npixels
        self.npix = hp.nside2npix(nside=nside)
        colors = ['lightgrey', 'lightgreen', 'red', 'lightblue']
        colors = ['gainsboro', 'lightgreen', 'red', 'blue']
        self.minx = 0
        self.maxx = len(colors)
        self.norm = plt.cm.colors.Normalize(self.minx, self.maxx)
        self.cmap = mpl.colors.ListedColormap(colors)
        self.cmap.set_under('w')
        self.outDir = outDir

    def visu(self, tab, vardisp='color', healpixId='healpixID', title='', inum=1, save=False):
        """
        Method to display OS (Mollview)

        Parameters
        --------------
        tab: array
          array of data
        vardisp: str, opt
           variable to display (default: color)
        healpixId: str, opt
          colname for healpixId (default: healpixID)
        title: str, opt
          title for the plot (default: )
        inum: int, opt
          tag for outputfile (default: -1)

        """
        fig = plt.figure()
        hpxmap = np.zeros(self.npix, dtype=np.float)
        hpxmap += -1
        healpixIDs = np.unique(tab[healpixId])
        hpxmap[tab['healpixID'].astype(int)] = tab[vardisp]

        hp.mollview(hpxmap, fig=fig, min=self.minx, max=self.maxx, cmap=self.cmap,
                    title=title, nest=True, norm=self.norm)
        hp.graticule(dpar=10., dmer=30.)

        if save:
            figName = '{}/obs_{}'.format(self.outDir, inum)
            plt.savefig(figName)

        # fig.colorbar(self.cmap)
        plt.close(fig)
        # plt.show()


class PlotNights(PlotOS):
    """
    class to plot all the nights for OS - inherit from plotOS

    Parameters
    --------------
    data: array
      data to display
    nside: int, opt
      nside for healpix (default: 128)
    outDir: str, opt
      output directory (default: plots_cadence)
    nproc: int, opt
      number pf procs for multiprocessing (default: 8)


    """

    def __init__(self, data, nside=128, outDir='plots_cadence', nproc=8):
        super().__init__(nside, outDir)

        if 'night' not in data.columns:
            data['night'] = data['time']-data['time'].min()+1
            data['night'] = data['night'].astype(int)

        # add color col
        data['color'] = 1.5
        ida = data['band'] == 'ztfr'
        data.loc[ida, 'color'] = 2.5
        ida = data['band'] == 'ztfi'
        data.loc[ida, 'color'] = 3.5

        self.data = data
        self.outDir = outDir
        checkDir(outDir)
        self.nproc = nproc

    def __call__(self, run_type='per_night'):

        night_max = self.data['night'].max()
        nights = range(1, night_max+1)
        # var = 'night'

        """
        if run_type == 'per_exposure':
            nights = self.data['time'].unique()
            var = 'time'
        """

        params = {}
        params['run_type'] = run_type
        multiproc(nights, params, self.processNight, self.nproc)

        """
        for night in nights:

            idx = self.data['night'] <= night
            # obs_night = self.transform(self.data[idx])
            obs_night = self.data[idx]
            if run_type == 'per_night':
                self.plot_night(obs_night, night)

            if run_type == 'per_exposure':
                times = obs_night['time'].unique()
                for ti in times:
                    io = obs_night['time'] <= ti
                    sel = obs_night[io]
                    self.plot_night(sel, 'time')
        """

    def plot_night(self, obs, night, var='night'):

        if len(obs) == 0:
            obs = pd.DataFrame(
                [(1, 0)], columns=['healpixID', 'color'])
        """
        else:
            ido = obs[var] < night
            if len(obs[ido]) > 0:
                obs.loc[ido, 'color'] = 0.5
        """
        tit = 'night {}'.format(night)
        if var != 'night':
            tit += ' - JD {}'.format(np.round(mjd, 2))
        self.visu(obs.to_records(index=False),
                  title=tit, inum=night, save=True)

    def transform(self, data):

        df = pd.DataFrame()

        # get healpixIDs
        hIDs = ','.join(data['healpixID'].to_list())
        hIDs = set(hIDs.split(','))
        if len(hIDs) == 0:
            return df
        print(hIDs)
        if 'None' in hIDs:
            hIDs.remove('None')

        for hID in hIDs:
            idx = data['healpixID'].str.contains(hID)
            dd = pd.DataFrame(data[idx])
            dd['healpixID'] = hID
            df = pd.concat((df, dd))

        return df

    def processNight(self, nights, params, j=0, output_q=None):

        run_type = params['run_type']
        for night in nights:
            idx = self.data['night'] == night
            obs_tonight = pd.DataFrame(self.data.loc[idx])
            obs_tonight = pd.DataFrame(obs_tonight[['healpixID', 'color']])
            obs_before = self.obs_night_before(
                night, obs_tonight['healpixID'].to_list())

            # obs_night = self.transform(self.data[idx])

            obs_night = pd.concat((obs_tonight, obs_before))

            if run_type == 'per_night':
                self.plot_night(obs_night, night)

            if run_type == 'per_exposure':
                times = obs_night['time'].unique()
                for ti in times:
                    io = obs_night['time'] <= ti
                    sel = obs_night[io]
                    self.plot_night(sel, 'time')

        if output_q is not None:
            return output_q.put({j: 1})
        else:
            return 1

    def obs_night_before(self, night, pixelList):

        ido = self.data['night'] < night
        obs_before = pd.DataFrame(self.data[ido])
        if len(obs_before) > 0:
            obs_before = pd.DataFrame(
                obs_before['healpixID'].unique(), columns=['healpixID'])
            obs_before['color'] = 0.5
            idx = obs_before['healpixID'].isin(pixelList)
            obs_before = obs_before[~idx]

        return obs_before
