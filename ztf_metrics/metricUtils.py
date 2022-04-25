from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery
import healpy as hp
import pandas as pd


def dustMap(nside=64):
    """
    Method to get the dust map

    Parameters
    ---------------
    nside: int, opt
      healpix nside parameters

    Returns
    -----------
    dustmap: pandas df
    with the following columns: healpixID, RA, Dec, ebvofMW

    """

    # get the total number of pixels
    npixels = hp.nside2npix(nside)

    # pixels = hp.get_all_neighbours(nside, 0.0, 0.0, nest=True, lonlat=Tru)

    # get the (RA, Dec) of the pixel centers
    vec = hp.pix2ang(nside, range(npixels), nest=True, lonlat=True)

    dustmap = pd.DataFrame(range(npixels), columns=['healpixID'])
    dustmap['RA'] = vec[0]
    dustmap['Dec'] = vec[1]

    coords = SkyCoord(vec[0], vec[1], unit='deg')
    try:
        sfd = SFDQuery()
    except Exception as err:
        from dustmaps.config import config
        config['data_dir'] = 'dustmaps'
        import dustmaps.sfd
        dustmaps.sfd.fetch()
    sfd = SFDQuery()
    ebvofMW = sfd(coords)

    dustmap['ebvofMW'] = ebvofMW
    dustmap['npixels'] = npixels

    return dustmap


def coaddNight(data, cols=['night', 'band', 'healpixID']):
    """
    Function to perform coaddition per band and per night

    Parameters
    ---------------
    data: pandas df
      data to process
    cols: list(str), opt
      list of columns for the groupby (default: ['night','band','healpixID'])

    Returns
    ----------
    pandas df of the mean from the groupby by cols

    """
    res = data.groupby(cols).median().reset_index()
    return res


def seasons(df, gap=60, mjdCol='time'):
    """
    Function that determines the different seasons of an observation related to a pixel.

    Parameters
    --------------
    df: pandas df
        data to process with observations.
    gap: int, opt
      gap between two seasons (default : 60 days)
    mjdCol: str, opt
      name of the mjd column (default: time)

    Returns
    --------------
    Starting pandas with a 'season' column.
    """

    df = df.sort_values(mjdCol, ascending=True)
    df = df.reset_index()
    df_time = df[mjdCol].diff()
    index = list((df_time[df_time > gap].index)-1)
    index.insert(0, df.index.min())
    index.insert(len(index), df.index.max())
    df['season'] = 0

    for i in range(0, len(index)-1):
        df.loc[index[i]: index[i+1], 'season'] = i+1
    return df


def addNight(data, mjdCol='time'):
    """
    function to add a night column to df

    Parameters
    --------------
    data: pandas df
      data to process
    mjdCol: str, opt
      name of the MJD col

    Returns
    ----------
    original pandas df + night column
    """
    data['night'] = data[mjdCol]-data[mjdCol].min()
    data['night'] = data['night'].astype(int)+1

    return data
