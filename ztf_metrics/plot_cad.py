import matplotlib.pyplot as plt
import healpy as hp
import matplotlib as mpl
import numpy as np
from ztf_pipeutils.ztf_util import checkDir, multiproc
import pandas as pd
from matplotlib import cm


class PlotCad:

    def __init__(self, tab=None, color=['green', 'red', 'black'], ls= [':', '-.', '--'],cad_min=0, cad_max=6):

        self.tab = tab
        self.n = np.arange(cad_min,cad_max,1)
        self.color = color
        self.ls = ls
        
    def hist_cad(self, all_band=False):
        ratio = [] 
        cad_ = []
        band_ = []
        for cad in self.n:
            for band in ['cad_ztfg', 'cad_ztfr', 'cad_ztfi']:
                mask = (self.tab[band]>=cad) & (self.tab[band]<cad+1)
                tab_cad = self.tab[mask]
                ratio.append((tab_cad['healpixID'].count()/self.tab['healpixID'].count()))
                cad_.append(cad+0.5)
                band_.append(band)
        
        
        df = pd.DataFrame({'cad': cad_, 'ratio': ratio, 'band': band_})

        df_g = df[df['band']=='cad_ztfg']
        df_r = df[df['band']=='cad_ztfr']
        df_i = df[df['band']=='cad_ztfi']

        add = [-0.25, 0.0, 0.25]
        labels = ['ac g band', 'ac r band', 'ac i band']

        fig, ax1 = plt.subplots(1, 1, figsize=(8,6))
        for i, df_ in enumerate ([df_g, df_r, df_i]):
            df_['accumulate_ratio'] = np.add.accumulate(df_['ratio'])
            bar = ax1.bar(x = df_['cad']+add[i], height = df_['ratio'], width=0.25, alpha=0.5, color=self.color[i],
                         label='{}'.format(df_['band'].unique()[0]))
    
            p, = ax1.plot(df_['cad']+add[i], df_['accumulate_ratio'], ls=self.ls[i], marker='o', color=self.color[i], label=labels[i])

        leg = ax1.legend(title='Accumulate ratio ("ac")\n& cadence :', loc='best')

        ax1.set_xlabel('Cadence')
        ax1.set_ylabel(r'$N_{SN}/N_{total}$')
        
        if all_band:
            n = np.arange(0,6,1)
            ratio = [] 
            cad_ = []
            band_ = []
            for cad in n:
                mask = (self.tab['cad_all']>=cad) & (self.tab['cad_all']<cad+1)
                tab_cad = self.tab[mask]
                ratio.append((tab_cad['healpixID'].count()/self.tab['healpixID'].count()))
                cad_.append(cad+0.5)
                band_.append('cad_all')
        
            df = pd.DataFrame({'cad': cad_, 'ratio': ratio, 'band': band_})

            fig2, ax1 = plt.subplots(1, 1, figsize=(8,6))
            df['accumulate_ratio'] = np.add.accumulate(df['ratio'])
            bar = ax1.bar(x = df['cad'], height = df['ratio'], width=0.25, alpha=0.5,
                 label='all')
    
            p, = ax1.plot(df['cad'], df['accumulate_ratio'], ls=':', marker='o', label='accumulate ratio')

            leg = ax1.legend(title='Accumulate ratio ("ac")\n& cadence :', loc='best')

            ax1.set_xlabel('Cadence')
            ax1.set_ylabel(r'$N_{SN}/N_{total}$')
            plt.show()
        
            return fig2
        
        return fig
        
    def plot_cad(self):
        new_tab = self.tab[self.tab['cad_ztfi']!=0]
        
        y = new_tab['cad_ztfi']
        x_g = new_tab['cad_ztfg']
        x_r = new_tab['cad_ztfr']
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,6))
        ax1.scatter(x_g, y, color='green')
        ax2.scatter(x_r, y, color='red')
        
        ax1.set_xlabel('Cadence band g')
        ax2.set_xlabel('Cadence band r')
        
        ax1.set_ylabel('Cadence band i')
        ax2.set_ylabel('Cadence band i')
        
        fig.tight_layout()
        
        return fig
        
       
        
        
        
        
        