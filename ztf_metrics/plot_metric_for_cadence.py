import matplotlib.pyplot as plt
import healpy as hp
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import cm
import seaborn as sns
import os


class PlotCAD:

    def __init__(self, nside=128, outDir='plot_cadence'):

        # get npixels
        self.npix = hp.nside2npix(nside=nside)
        self.colors = ['white', 'green', 'red', 'orange',
                       'turquoise', 'dimgray', 'dimgray', 'dimgray']
        self.minx = 0
        self.maxx = len(self.colors)
        self.norm = plt.cm.colors.Normalize(self.minx, self.maxx)
        self.cmap = mpl.colors.ListedColormap(self.colors)
        self.outDir = outDir
        self.pixArea = hp.nside2pixarea(nside, degrees=True)

    def visu(self, tab, vardisp, healpixId='healpixID', title='', inum=1, save=False, cbar=True):
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

        self.sel = tab
        print('moyenne =', np.mean(self.sel[vardisp]))
        print('medianne =', np.median(self.sel[vardisp]))

        hpxmap = np.zeros(self.npix, dtype=np.float)
        hpxmap += -1
        healpixIDs = np.unique(self.sel[healpixId])
        hpxmap[self.sel['healpixID'].astype(int)] = self.sel[vardisp]

        ik = hpxmap > 1.
        area = self.pixArea*len(hpxmap[ik])

        title += ''

        hp.mollview(hpxmap, fig=fig, min=self.minx, max=self.maxx, cmap=self.cmap,
                    title=title, nest=True, norm=self.norm, cbar=cbar)
        hp.graticule(dpar=10., dmer=30.)

        # define a new colorbar (no access to mollview one)
        cax = fig.add_axes(
            [ax.get_position().x0+0.2, ax.get_position().y0+0.05, ax.get_position().width-0.4, 0.05])
        bounds = [0, 1, 2, 3, 4, 5, np.max(self.sel[vardisp])]
        cb2 = mpl.colorbar.ColorbarBase(cax, cmap=self.cmap,
                                        norm=self.norm,
                                        boundaries=bounds,
                                        ticks=bounds,
                                        orientation='horizontal')
        cb2.set_label('test')

        if not os.path.exists(self.outDir):
            os.makedirs(self.outDir)
            figName = '{}/_{}'.format(self.outDir, title)
            fig.savefig(figName)
        else:
            figName = '{}/_{}'.format(self.outDir, title)
            fig.savefig(figName)

        return fig
