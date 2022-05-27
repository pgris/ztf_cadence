import matplotlib.pyplot as plt
import healpy as hp
import matplotlib as mpl
import numpy as np
from ztf_pipeutils.ztf_util import checkDir, multiproc
import pandas as pd
from matplotlib import cm
import seaborn as sns


class PlotOS:
    """
    class to plot OS pointings vs time
    Parameters
    --------------
    nside: int, opt
      nside for healpix (default: 128)
    outDir: str,opt
      output directory (default: plot_cadence)
    hold: bool, opt
      to make multi var on the same plot (default: False)
    """

    def __init__(self, nside=128, outDir='plot_cadence', hold=False, nround=2):

        # get npixels
        self.npix = hp.nside2npix(nside=nside)
        self.cmap = mpl.cm.get_cmap('jet').copy()
        self.cmap.set_extremes(under='w', over='gray')
        self.outDir = outDir
        self.pixArea = hp.nside2pixarea(nside, degrees=True)
        self.hold = hold
        self.nround = nround

    def visu(self, tab, vardisp='color', healpixId='healpixID', title='', inum=1, save=False, min_tab=0, max_tab=100,
             cbar=True):
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

        self.cbar = cbar
        idx = tab['healpixID'] == 'None'
        self.sel = tab[~idx]
        self.min_ = min_tab
        self.max_ = max_tab
        norm = plt.cm.colors.Normalize(self.min_, self.max_)

        hpxmap = np.zeros(self.npix, dtype=np.float)
        hpxmap += -1
        healpixIDs = np.unique(self.sel[healpixId])
        hpxmap[self.sel['healpixID'].astype(int)] = self.sel[vardisp]

        mk_sup = self.sel[vardisp] < max_tab
        mk_inf = self.sel[vardisp] > min_tab
        mk_born = mk_sup & mk_inf

        tab_masked = self.sel[mk_born]

        title += ' / median = {} \n moy = {}'.format(
            round(tab_masked[vardisp].median(), self.nround), round(tab_masked[vardisp].mean(), self.nround))

        ik = hpxmap > 1.
        area = self.pixArea*len(hpxmap[ik])

        hp.mollview(hpxmap, min=self.min_, max=self.max_, cmap=self.cmap,
                    title=title, nest=True, norm=norm, cbar=self.cbar, hold=self.hold)
        hp.graticule(dpar=10., dmer=30.)
