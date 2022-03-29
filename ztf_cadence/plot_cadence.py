import matplotlib.pyplot as plt
import healpy as hp
import matplotlib as mpl
import numpy as np
from ztf_pipeutils.ztf_util import checkDir


class PlotOS:
    """
    class to plot OS pointings vs time

    Parameters
    --------------
    nside: int, opt
      nside for healpix (default: 128)


    """

    def __init__(self, nside=128, outDir='plot_cadence'):

        # get npixels
        self.npix = hp.nside2npix(nside=nside)
        colors = ['grey', 'lightgreen', 'red', 'lightblue']
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
        hpxmap = np.zeros(self.npix, dtype=np.float)
        hpxmap += -1
        healpixIDs = np.unique(tab[healpixId])
        hpxmap[tab['healpixID'].astype(int)] = tab[vardisp]

        hp.mollview(hpxmap, min=self.minx, max=self.maxx, cmap=self.cmap,
                    title=title, nest=True, norm=self.norm)
        hp.graticule(dpar=10., dmer=30., verbose=False)

        if save:
            figName = '{}/obs_{}'.format(self.outDir, inum)
            plt.savefig(figName)

            # plt.show()


class PlotNights(PlotOS):
    """
    class to plot all the nights for OS - inherit from plotOS

    Parameters
    --------------
    data: array
      data to display
    nside: int, opt
      nside for healpix

    """

    def __init__(self, data, nside=128, outDir='plots_cadence'):
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

    def __call__(self, run_type='per_night'):

        night_max = self.data['night'].max()
        nights = range(1, night_max+1)
        #var = 'night'

        """
        if run_type == 'per_exposure':
            nights = self.data['time'].unique()
            var = 'time'
        """

        for night in nights:

            idx = self.data['night'] <= night
            obs_night = self.data[idx]
            if run_type == 'per_night':
                self.plot_night(obs_night, night)

            if run_type == 'per_exposure':
                times = obs_night['time'].unique()
                for ti in times:
                    io = obs_night['time'] <= ti
                    sel = obs_night[io]
                    self.plot_night(sel, 'time')

    def plot_night(self, obs, night, var='night'):

        if len(obs) == 0:
            obs = pd.DataFrame(
                [(1, 0)], columns=['healpixID', 'color'])
        else:
            ido = obs[var] < night
            if len(obs[ido]) > 0:
                obs.loc[ido, 'color'] = 0.5

        tit = 'night {}'.format(night)
        if var != 'night':
            tit += ' - JD {}'.format(np.round(mjd, 2))
        self.visu(obs.to_records(index=False),
                  title=tit, inum=night, save=True)
