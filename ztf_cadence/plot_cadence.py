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
        self.colors = ['gainsboro', 'green', 'red', 'magenta']
        self.minx = 0
        self.maxx = len(self.colors)
        self.norm = plt.cm.colors.Normalize(self.minx, self.maxx)
        self.cmap = mpl.colors.ListedColormap(self.colors)
        self.cmap.set_under('w')
        self.outDir = outDir
        self.pixArea = hp.nside2pixarea(nside, degrees=True)

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
        #fig, ax = plt.subplots()
        ax = fig.add_axes([0, 0, 1, 1])
        hpxmap = np.zeros(self.npix, dtype=np.float)
        hpxmap += -1
        healpixIDs = np.unique(tab[healpixId])
        hpxmap[tab['healpixID'].astype(int)] = tab[vardisp]

        ik = hpxmap > 1.
        area = self.pixArea*len(hpxmap[ik])

        title += ' - area {} deg$^2$'.format(int(area))
        hp.mollview(hpxmap, fig=fig, min=self.minx, max=self.maxx, cmap=self.cmap,
                    title=title, nest=True, norm=self.norm, cbar=False)
        hp.graticule(dpar=10., dmer=30.)

        # define a new colorbar (no access to mollview one)
        cax = fig.add_axes(
            [ax.get_position().x0+0.2, ax.get_position().y0+0.05, ax.get_position().width-0.4, 0.05])
        bounds = range(0, 5)
        cb2 = mpl.colorbar.ColorbarBase(cax, cmap=self.cmap,
                                        norm=self.norm,
                                        boundaries=bounds,
                                        ticks=bounds,
                                        orientation='horizontal',
                                        ticklocation='bottom')
        cb2.ax.get_yaxis().set_ticks([])
        cb2.ax.get_xaxis().set_ticks([])
        for j, lab in enumerate(['previous', 'g', 'r', 'i']):
            k = 0.5
            if j == 0:
                k = 0.3
            cb2.ax.text(k+j, 0.80, lab)  # ha='center', va='center')

        # hide axis
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.axis('off')

        if save:
            figName = '{}/obs_{}'.format(self.outDir, inum)
            plt.savefig(figName)
        # plt.show()
        plt.close(fig)


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

    def processNight(self, nights, params, j=0, output_q=None):

        run_type = params['run_type']
        for night in nights:
            obs_night = self.obs_tonight(night)
            pixels = []
            if len(obs_night) > 0:
                pixels = obs_night['healpixID'].to_list()
            obs_before = self.obs_night_before(night, pixels)

            # obs_night = self.transform(self.data[idx])

            obs_night = pd.concat((obs_night, obs_before))

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

    def obs_tonight(self, night):

        idx = self.data['night'] == night
        obs_night = self.data[idx]
        obs = pd.DataFrame()
        for b in obs_night['color'].unique():
            io = obs_night['color'] == b
            sel = obs_night[io]
            dd = pd.DataFrame(self.get_hID(sel), columns=['healpixID'])
            dd['color'] = b
            obs = pd.concat((obs, dd))

        return obs

    def obs_night_before(self, night, pixelList):

        ido = self.data['night'] < night
        obs_before = pd.DataFrame(self.data[ido])
        if len(obs_before) > 0:
            obs_before = pd.DataFrame(
                self.get_hID(obs_before), columns=['healpixID'])
            obs_before['color'] = 0.5
            idx = obs_before['healpixID'].isin(pixelList)
            obs_before = obs_before[~idx]

        return obs_before

    def get_hID(self, data):

        hIDs = ','.join(data['healpixID'].to_list())
        hIDs = set(hIDs.split(','))
        if 'None' in hIDs:
            hIDs.remove('None')

        return hIDs
