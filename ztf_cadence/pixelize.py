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

    """

    def __init__(self, nside=128, width_ra=0.85, width_dec=0.85):

        self.nside = nside
        self.width_ra = width_ra
        self.width_dec = width_dec

    def __call__(self, df):

        dftot = pd.DataFrame()

        for i, dd in df.iterrows():
            ppix = self.pixels(dd['ra'], dd['dec'])
            if len(ppix) > 0:
                dfa = pd.DataFrame(ppix, columns=['healpixID'])
                dfa['field'] = dd['field']
                dfa['field'] = dfa['field'].astype(int)
                dfa = dfa.merge(pd.DataFrame([dd]), left_on=[
                    'field'], right_on=['field'])
                dfa = dfa.drop_duplicates(subset=['healpixID'])

                dftot = pd.concat((dftot, dfa))

        return dftot

    def pixels(self, ra_center, dec_center):
        """"
        Method to grab the pixels corresponding to a polygon
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

        ra_min = ra_center-self.width_ra/2.
        dec_min = dec_center-self.width_dec/2.

        ra_max = ra_min+self.width_ra
        dec_max = dec_min+self.width_dec

        ras = [ra_min, ra_min, ra_max, ra_max]
        decs = [dec_min, dec_max, dec_max, dec_min]

        rr = pd.DataFrame(ras, columns=['ra'])
        rr['dec'] = decs
        arr = hp.ang2vec(rr['ra'], rr['dec'], lonlat=True)

        return hp.query_polygon(self.nside, arr, nest=True)
