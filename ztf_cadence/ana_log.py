import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
import matplotlib as mpl
import healpy as hp


def align(df):

    res = df.copy()
    col1_ra = [3, 0, 19, 16, 35, 32, 51, 48]
    col2_ra = [2, 1, 18, 17, 34, 33, 50, 49]

    idx1 = res['rcid'].isin(col1_ra)

    ra_mean = res[idx1]['ra'].mean()

    #width_ra = ra2_mean-ra1_mean
    width_ra = 0.85
    res.loc[idx1, 'ra'] = ra_mean

    for i in range(1, 3):
        bb = 4*i
        delta = 2*i*width_ra
        reslist = list(map(lambda x: x + bb, col1_ra))
        print('rr', reslist)
        idx = res['rcid'].isin(reslist)
        res.loc[idx, 'ra'] = ra_mean-delta*width_ra

    return res


def pixelize(df, outName, nside=128):
    width_ra = 0.89
    width_dec = 0.89

    dftot = pd.DataFrame()
    for i, dd in df.iterrows():
        ppix = pixels(dd['ra'], dd['dec'], width_ra, width_dec, nside)
        if len(ppix) > 0:
            dfa = pd.DataFrame(ppix, columns=['healpixID'])
            dfa['field'] = dd['field']
            dfa['field'] = dfa['field'].astype(int)
            dfa = dfa.merge(pd.DataFrame([dd]), left_on=[
                'field'], right_on=['field'])
            dfa = dfa.drop_duplicates(subset=['healpixID'])

        dftot = pd.concat((dftot, dfa))
        if (i % 10000) == 0:
            print(dftot)
            print('processed', i)
        if i > 100000:
            break

    dftot.to_hdf(outName, key='os')


def plotViewIndiv(nside, tab, xval, minx, maxx, healpixId='healpixID', title='', inum=1):
    """
    Mollweid view

    Parameters
    --------------

    nside: int
     healpix nside parameter
    tab: array
     array of values
    xval: str
     name of the variable to display
    legx: str
     legend to display
    unitx: str
     unit of the legend
    minx: float
     min value for display
    maxx: float
     max value for display
    band: str
     band considered
    dbName: str
     name of the cadence file
    saveFig: bool, opt
     To save figure or not (default: False)
    seasons: int
     list of seasons to display (-1 = all)
    type: str, opt
     type of display: mollview (default) or cartview
    fieldzoom: bool, opt
     to make a zoom centered on (fieldzoom['RA'],fieldzoom['Dec'])
     (default: None)
    healpixId: str, opt
     string to identify healpixId (default: healpixID)

    Returns
    ---------
    None

    """

    # fig, ax = plt.subplots(figsize=(8, 6))
    # print(xval,tab[xval])

    npix = hp.nside2npix(nside=nside)
    norm = plt.cm.colors.Normalize(minx, maxx)
    cmap = plt.cm.jet
    cmap = mpl.colors.ListedColormap(['lightgreen', 'red', 'lightblue'])

    # bounds = [-1.0, -0.5, 0.0, 0.5, 1.0]
    # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    # cmap = plt.cm.get_cmap('hsv', 3)
    cmap.set_under('w')
    # cmap.set_bad('w')

    hpxmap = np.zeros(npix, dtype=np.float)
    hpxmap += -1
    healpixIDs = np.unique(tab[healpixId])
    hpxmap[tab['healpixID'].astype(int)] = tab[xval]

    # plt.axes(ax)
    # plt.sca(ax)

    hp.mollview(hpxmap, min=minx, max=maxx, cmap=cmap,
                title=title, nest=True, norm=norm)
    hp.graticule(dpar=10., dmer=30., verbose=False)

    figName = 'figs_night/obs_{}'.format(inum)
    plt.savefig(figName)
    plt.show()


def select(df, colName, refVal):

    idx = np.abs(df[colName]-refVal) < 1.e-5
    return df[idx]


def pixels(ra_center, dec_center, width_ra, width_dec, nside=64):

    ra_min = ra_center-width_ra/2.
    dec_min = dec_center-width_dec/2.

    ra_max = ra_min+width_ra
    dec_max = dec_min+width_dec

    # config = [(ra_min, dec_min), (ra_min, dec_max),
    #          (ra_max, dec_max), (ra_max, dec_min)]
    ras = [ra_min, ra_min, ra_max, ra_max]
    decs = [dec_min, dec_max, dec_max, dec_min]

    rr = pd.DataFrame(ras, columns=['ra'])
    rr['dec'] = decs
    arr = hp.ang2vec(rr['ra'], rr['dec'], lonlat=True)

    return hp.query_polygon(nside, arr, nest=True)


def get_width(data, colName):

    selb = data.sort_values(by=[colName])

    min_val = data[colName].min()
    selb['diff'] = selb[colName]-min_val
    idx = np.abs(selb['diff']) > 0.2
    selc = selb[idx].reset_index()

    return selc.loc[0, 'diff'].mean()


def plot_cadence(night, df, nside=128):

    cols = dict(zip(['ztfg', 'ztfr', 'ztfi'], [0.5, 1.5, 2.5]))
    rpix = []

    selt = select(df, 'night', night)
    times = selt['time'].unique()
    dfpix = pd.DataFrame()
    for tti in times:
        sel = select(selt, 'time', tti)
        ra_center = sel['ra'].mean()
        dec_center = sel['dec'].mean()
        width_ra = np.sqrt(46.72)
        width_dec = width_ra
        pixlist = pixels(ra_center, dec_center,
                         width_ra, width_dec, nside=nside)
        band = sel['band'].unique()[0]
        dd = pd.DataFrame(pixlist, columns=['healpixID'])
        dd['color'] = cols[band]
        dfpix = pd.concat((dfpix, dd))
        """
        for pix in pixlist:
            rpix.append((pix, cols[band]))
        """
    # rr = np.rec.fromrecords(rpix, names=['healpixID', 'color'])

    if len(dfpix) > 0:
        rr = dfpix.to_records(index=False)
    else:
        rr = np.rec.fromrecords([(1, 0)],
                                names=['healpixID', 'color'])
    plotViewIndiv(nside, rr, 'color', 0, 3,
                  healpixId='healpixID', title='night {}'.format(night), inum=night)


# read the log file
outName = 'OS.hdf5'

fi = '../simu/data/2018_all_logs_from_dr1_rcid_zp_from_masci.csv'
df = pd.read_csv(fi)

"""
print(df.columns)
df = df.round({'time': 7})
time_ref = 2458197.6238542
il = df['time'] == time_ref
df = df[il]

dfb = align(df)

fig, ax = plt.subplots()
ax.plot(df['ra'], df['dec'], 'ko', mfc='None')
ax.plot(dfb['ra'], dfb['dec'], 'r*')

for rcid in range(64):
    idx = df['rcid'] == rcid
    sel = df[idx]
    #ax.plot(sel['ra'], sel['dec'], 'r*')
    wdth = 0.1
    ax.text(sel['ra']+wdth, sel['dec']+wdth, '{}'.format(int(rcid)))
"""
# plt.show()


nside = 128

#pixelize(df, outName, nside=nside)

df = pd.read_hdf(outName)
# add night col
df['night'] = df['time']-df['time'].min()+1
df['night'] = df['night'].astype(int)
print(df)
nights = df['night'].unique()
# add color col
df['color'] = 0.5
ida = df['band'] == 'ztfr'
df.loc[ida, 'color'] = 1.5
ida = df['band'] == 'ztfi'
df.loc[ida, 'color'] = 2.5
pixArea = hp.nside2pixarea(nside, degrees=True)
for night in range(1, 365):
    idx = df['night'] == night
    obs_night = df[idx]
    if len(obs_night) == 0:
        obs_night = pd.DataFrame([(1, 0)], columns=['healpixID', 'color'])
    else:
        times = obs_night['time'].unique()
        for ti in times:
            io = obs_night['time'] == ti
            sel = obs_night[io]
            rcids = sel['rcid'].unique()
            print('area', ti, len(sel)*pixArea, rcids,
                  len(rcids), len(sel), pixArea)
            """
            for rcid in rcids:
                print('quadrant', rcid, len(sel[sel['rcid'] == rcid]))
            """
    plotViewIndiv(nside, obs_night.to_records(index=False), 'color', 0, 3,
                  healpixId='healpixID', title='night {}'.format(night), inum=night)


df = df.round({'time': 7})
times = df['time'].unique()
time_min = df['time'].min()
df['night'] = df['time']-time_min+1
df['night'] = df['night'].astype(int)
nights = df['night'].unique()
print('ahhh', df.columns)
df = df.sort_values(by=['night'])
rt = []
for night in range(1, 365):
    plot_cadence(night, df)
    if night > 10:
        break
    # plt.show()
"""
field = 509.
sel = select(df, 'time', time_ref)
# sel = select(df, 'field', field)
print(len(sel), sel['band'].unique())

ax.plot(sel['ra'], sel['dec'], 'k.')

diff_ra = get_width(sel, 'ra')
diff_dec = get_width(sel, 'dec')
print(diff_ra, diff_dec)
diff_ra /= 2.
diff_dec /= 2.
print(diff_ra, diff_dec)

ra_center = sel['ra'].mean()
dec_center = sel['dec'].mean()
ax.plot([ra_center], [dec_center], 'ks')
width_ra = 75.32 - 67.917
width_dec = 15.55-8.065
width_ra = np.sqrt(46.72)
width_dec = width_ra
config, pixlist = pixels(ra_center, dec_center,
                         width_ra, width_dec, nside=nside)
print('aaera', width_ra, width_dec, width_ra*width_dec)
rpix = []
for pix in pixlist:
    rpix.append((pix, 0.5))
"""
"""
for i, row in sel.iterrows():
    ra = row['ra']-np.abs(diff_ra)
    dec = row['dec']-np.abs(diff_dec)
    print(row['ra'], row['dec'], ra, dec)
    rect = patches.Rectangle((ra, dec), 2.*diff_ra, 2.*diff_dec,
                             linewidth=1, edgecolor='r', facecolor='none')
    # Add the patch to the Axes
    ax.add_patch(rect)
    print('area', rect.get_height()*rect.get_width())
    config, pixlist = pixels(
        row['ra'], row['dec'], 2.*diff_ra, 2.*diff_dec, nside=nside)
    cc = pd.DataFrame(config, columns=['ra', 'dec'])
    ax.plot(cc['ra'], cc['dec'], 'r*')
    print(pixlist)
    for pix in pixlist:
        rpix.append((pix, 0.5))
    # break
"""
# rr = np.rec.fromrecords(rpix, names=['healpixID', 'color'])
# plotViewIndiv(nside, rr, 'color', 0.01, 1., healpixId='healpixID')

"""
for i, row in sel.iterrows():
    print(row[['ra', 'dec']])
    print(pixels(row['ra'], row['dec'], 2.*diff_ra, 2.*diff_dec))
"""
# plt.show()
