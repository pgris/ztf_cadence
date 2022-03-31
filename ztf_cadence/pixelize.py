import healpy as hp
import pandas as pd


class Pixelize_sky:
    """
    class to project observations on the pixelized sky (with healpix)

    Parameters
    --------------
    nside: int, opt
      nside parameter for healpix (default: 128)
    width_ra: float, opt
       ra width of a quadrant (default: 0.85)
    width_dec: float, opt
       dec width of a quadrant (default: 0.85)
    raCol: str, opt
      name of the RA col for pixelizing (default: ra)
    decCol: str, opt
      name of the Dec col for pixelizing (default: dec)

    """

    def __init__(self, nside=128, width_ra=0.85, width_dec=0.85, raCol='ra', decCol='dec'):

        self.nside = nside
        self.width_ra = width_ra
        self.width_dec = width_dec
        self.raCol = raCol
        self.decCol = decCol

    def __call__(self, df):

        dftot = pd.DataFrame()
        import time
        for i, dd in df.iterrows():
            time_ref = time.time()
            ppix = pixels(self.nside, dd[self.raCol], dd[self.decCol],
                          self.width_ra, self.width_dec)
            """
            if len(ppix) > 0:
                 dfa = pd.DataFrame(ppix, columns=['healpixID'])
                dfa['field'] = dd['field']
                dfa['field'] = dfa['field'].astype(int)
                dfa = dfa.merge(pd.DataFrame([dd]), left_on=[
                    'field'], right_on=['field'])
                #dfa = dfa.drop_duplicates(subset=['healpixID'])
            """
            dfa = pd.DataFrame([dd])
            if len(ppix) > 0:
                dfa['healpixID'] = ','.join(map(str, ppix))
            else:
                dfa['healpixID'] = ['None']
            #print('there man', dfa)
            dftot = pd.concat((dftot, dfa))
        return dftot


def pixels(nside, ra_center, dec_center, width_ra, width_dec):
    """"
    function to grab the pixels corresponding to a polygon
    with the center (ra_center, dec_center) and width_ra, width_dec as
    widths in (ra,dec)

    Parameters
    ---------------
    ra_center: float
       ra value of the center
    dec_center: float
      dec value of the center

    Returns
    ----------
    list of healpixIDs inside the polygon


    """

    ra_min = ra_center-width_ra/2.
    dec_min = dec_center-width_dec/2.

    ra_max = ra_min+width_ra
    dec_max = dec_min+width_dec

    ras = [ra_min, ra_min, ra_max, ra_max]
    decs = [dec_min, dec_max, dec_max, dec_min]

    rr = pd.DataFrame(ras, columns=['ra'])
    rr['dec'] = decs
    arr = hp.ang2vec(rr['ra'], rr['dec'], lonlat=True)

    return hp.query_polygon(nside, arr, nest=True)


class Pixel2Obs:

    def __init__(self, nside, data, width_ra, width_dec, raCol='ra', decCol='dec'):

        self.nside = nside
        self.data = data
        self.wr = 1.1*width_ra/2.
        self.wd = 1.1*width_dec/2.
        self.raCol = raCol
        self.decCol = decCol

        self.pixelize_it = Pixelize_sky(nside, width_ra=width_ra,
                                        width_dec=width_dec, raCol=raCol, decCol=decCol)

    def __call__(self, healpixID):

        ra_pixel, dec_pixel = hp.pix2ang(
            self.nside, healpixID, nest=True, lonlat=True)

        datasel = self.select(ra_pixel, dec_pixel)

        if datasel.empty:
            return pd.DataFrame()
        else:
            dd = self.pixelize_it(datasel)
            if not dd.empty:
                it = dd['healpixID'] == healpixID
                return dd[it]
            else:
                return pd.DataFrame()

    def select(self, ra_pixel, dec_pixel):

        rmin = ra_pixel-self.wr
        rmax = ra_pixel+self.wr
        dmin = dec_pixel-self.wd
        dmax = dec_pixel+self.wd
        idx = self.data[self.raCol] >= rmin
        idx &= self.data[self.raCol] < rmax
        idx &= self.data[self.decCol] >= dmin
        idx &= self.data[self.decCol] < dmax

        return self.data[idx]
