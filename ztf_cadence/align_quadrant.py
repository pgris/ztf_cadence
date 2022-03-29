import pandas as pd
import numpy as np
from sklearn.cluster import KMeans


class Cluster:

    def __init__(self, data, nclusters, RA_name='ra', Dec_name='dec'):
        """
        class to identify clusters of points in (RA,Dec)
        Parameters
        ---------------
        data: numpy record array
         data to process
        nclusters: int
         number of clusters to find
        RA_name: str, opt
         field name for the RA (default=fieldRA)
        Dec_name: str, opt
         field name for the Dec (default=fieldDec)
        """

        # grab necessary infos
        self.data = data
        self.RA_name = RA_name
        self.Dec_name = Dec_name

        # make the cluster of points
        self.points, self.clus, self.labels = self.makeClusters(nclusters)

        # analyse the clusters
        clusters, dfclus, dataclus = self.anaClusters(nclusters)

        # this is a summary of the clusters found
        self.clusters = clusters
        self.dfclusters = dfclus
        self.dataclus = dataclus

        """
        # this dataframe is a matching of initial data and clusters infos
        datadf = pd.DataFrame(np.copy(data))
        print(datadf.columns)
        print(dfclus[['RA', 'Dec']])
        datadf[self.RA_name] = datadf[self.RA_name].round(3)
        datadf[self.Dec_name] = datadf[self.Dec_name].round(3)
        dfclus['RA'] = dfclus['RA'].round(3)
        dfclus['Dec'] = dfclus['Dec'].round(3)
        self.dataclus = datadf.merge(
            dfclus, left_on=[self.RA_name, self.Dec_name], right_on=['RA', 'Dec'])
        """

    def makeClusters(self, nclusters):
        """
        Method to identify clusters
        It uses the KMeans algorithm from scipy
        Parameters
        ---------------
        nclusters: int
         number of clusters to find
        Returns
        -----------
        points: numpy array
          array of (RA,Dec) of the points
        y_km: numpy array
          index of the clusters
        kmeans.labels_: numpy array
          kmeans label
        """

        """
        r = []
        for (pixRA, pixDec) in self.data[[self.RA_name,self.Dec_name]]:
            r.append([pixRA, pixDec])
        points = np.array(r)
        """

        points = np.array(
            self.data[[self.RA_name, self.Dec_name]].tolist())
        # points = np.array(self.data[[self.RA_name, self.Dec_name]
        #                            ].to_records(index=False))
        # create kmeans object
        kmeans = KMeans(n_clusters=nclusters)
        # fit kmeans object to data
        kmeans.fit(points)

        # print location of clusters learned by kmeans object
        # print('cluster centers', kmeans.cluster_centers_)

        # save new clusters for chart
        y_km = kmeans.fit_predict(points)

        return points, y_km, kmeans.labels_

    def anaClusters(self, nclusters):
        """
        Method matching clusters to data
        Parameters
        ---------------
        nclusters: int
         number of clusters to consider
        Returns
        -----------
        env: numpy record array
          summary of cluster infos:
          clusid, fieldId, RA, Dec, width_RA, width_Dec,
          area, dbName, fieldName, Nvisits, Nvisits_all,
          Nvisits_u, Nvisits_g, Nvisits_r, Nvisits_i,
          Nvisits_z, Nvisits_y
        dfcluster: pandas df
          for each data point considered: RA,Dec,fieldName,clusId
        """

        rcluster = pd.DataFrame()
        dfcluster = pd.DataFrame()
        datacluster = pd.DataFrame()
        for io in range(nclusters):
            RA = self.points[self.clus == io, 0]
            Dec = self.points[self.clus == io, 1]
            dfclus = pd.DataFrame({'RA': RA, 'Dec': Dec})
            # ax.scatter(RA,Dec, s=10, c=color[io])
            indx = np.where(self.labels == io)[0]
            sel_obs = self.data[indx]

            min_RA = np.min(RA)
            max_RA = np.max(RA)
            min_Dec = np.min(Dec)
            max_Dec = np.max(Dec)
            mean_RA = np.mean(RA)
            mean_Dec = np.mean(Dec)

            datacluster_loc = pd.DataFrame(np.copy(sel_obs))
            datacluster_loc.loc[:, 'clusId'] = int(io)
            datacluster_loc.loc[:, 'RA'] = mean_RA
            datacluster_loc.loc[:, 'Dec'] = mean_Dec
            datacluster = pd.concat([datacluster, datacluster_loc], sort=False)

            dfclus.loc[:, 'clusId'] = int(io)
            dfcluster = pd.concat([dfcluster, dfclus], sort=False)

            rclus = pd.DataFrame(columns=['clusid'])
            rclus.loc[0] = int(io)
            rclus.loc[:, 'RA'] = mean_RA
            rclus.loc[:, 'Dec'] = mean_Dec
            rclus.loc[:, 'width_RA'] = max_RA-min_RA
            rclus.loc[:, 'width_Dec'] = max_Dec-min_Dec

            rcluster = pd.concat((rcluster, rclus))

        return rcluster, dfcluster, datacluster


class Align_Field:
    """
    class to compute quadrant center position for a set of observations

    Parameters
    ---------------
    data: array
      data to process
    delta_ccd: float, opt
       ccd widht (in deg) (default: 1.71375)
    raCol: str, opt
      RA col name (default: ra)
    decCol: str, opt
      Dec col name (default: dec)
    nclusters: int, opt
      number of clusters to make in the (ra,dec) plane (default: 64)
    """

    def __init__(self, data, delta_ccd=1.71375, raCol='ra', decCol='dec', nclusters=64, ccd_map={}):

        self.data = data
        self.delta_ccd = delta_ccd
        self.raCol = raCol
        self.decCol = decCol
        self.raquadCol = '{}_quad'.format(self.raCol)
        self.decquadCol = '{}_quad'.format(self.decCol)

        self.ra_mean = data[raCol].mean()
        self.dec_mean = data[decCol].mean()
        self.nclusters = nclusters
        self.ccd_map = ccd_map

        # get the quadrant pos
        #self.quadrant_pos = self.map_quadrant()

        res = self.map_ccd()
        quad_map = pd.DataFrame()
        for i, vv in res.iterrows():
            bb = self.map_quadrant(vv[self.raCol], vv[self.decCol],
                                   vv['ccd'], npix=1, deltapix=4)
            quad_map = pd.concat((quad_map, bb))

        quad_map = quad_map.rename(
            columns={self.raCol: self.raquadCol, self.decCol: self.decquadCol})
        for pp in ['ccd', 'rcid']:
            quad_map[pp] = quad_map[pp].astype(int)
            self.data[pp] = self.data[pp].astype(int)
        self.quadrants = quad_map

    def __call__(self):
        """
        Main method to add two columns to initial data: ra_quadrant and dec_quadrant

        """
        # merge data with quadrant positions
        # this is a method using clustering - not reliable
        """
        tab = pd.concat(
            (self.data[[self.raCol, self.decCol]], self.quadrant_pos))
        clusters = Cluster(tab[[self.raCol, self.decCol]].to_records(
            index=False), self.nclusters).dataclus

        clusters = pd.DataFrame(clusters)[[self.raCol, self.decCol, 'clusId']]
        cols = [self.raCol, self.decCol]
        mclus = self.quadrant_pos.merge(clusters, left_on=cols, right_on=cols)

        datab = self.data.merge(clusters, left_on=cols, right_on=cols)
        mclus = mclus.rename(
            columns={self.raCol: self.RAquad, self.decCol: self.Decquad})

        dataf = datab.merge(mclus, left_on=['clusId'], right_on=['clusId'])
        """

        rcids = self.data['rcid'].unique()
        idx = self.quadrants['rcid'].isin(rcids)
        dataf = self.data.merge(self.quadrants[idx], left_on=[
                                'rcid'], right_on=['rcid'])

        return dataf

    def map_quadrant_deprecated(self):
        """
        Method to redefine quadrant centers from the center of the field

        """

        r = []
        for i in range(-4, 4, 1):
            xv = self.ra_mean+(2*i+1)*self.delta_ccd/4.
            for j in range(-4, 4, 1):
                yv = self.dec_mean+(2*j+1)*self.delta_ccd/4.
                r.append((xv, yv))

        res = pd.DataFrame(r, columns=[self.raCol, self.decCol])

        return res

    def map_ccd(self, npix=2, deltapix=2):
        """
        Method to redefine ccd centers from the center of the field

        """

        r = []
        for i in range(-npix, npix, 1):
            xv = self.ra_mean+(2*i+1)*self.delta_ccd/deltapix
            for j in range(-npix, npix, 1):
                yv = self.dec_mean+(2*j+1)*self.delta_ccd/deltapix
                ccd_num = self.ccd_map[(i, j)]

                r.append((int(i), int(j), xv, yv, ccd_num))

        res = pd.DataFrame(
            r, columns=['i', 'j', self.raCol, self.decCol, 'ccd'])
        res['i'] = res['i'].astype(int)
        res['j'] = res['j'].astype(int)

        return res

    def map_quadrant(self, ra_mean, dec_mean, ccd, npix=1, deltapix=2):
        """
        Method to redefine quadrants from the center of the ccd

        """
        ref_quad = dict(
            zip([(-1, 0), (0, 0), (-1, -1), (0, -1)], [1, 0, 2, 3]))

        r = []
        for i in range(-npix, npix, 1):
            xv = ra_mean+(2*i+1)*self.delta_ccd/deltapix
            for j in range(-npix, npix, 1):
                yv = dec_mean+(2*j+1)*self.delta_ccd/deltapix
                quadnum = ref_quad[(i, j)]+4*(ccd-1)
                r.append((int(i), int(j), xv, yv, ccd, quadnum))

        res = pd.DataFrame(
            r, columns=['i', 'j', self.raCol, self.decCol, 'ccd', 'rcid'])
        res['i'] = res['i'].astype(int)
        res['j'] = res['j'].astype(int)

        return res


def map_ccd():
    """
    Function to map position of ccd (from the center of the FP) to ccd num

    """
    pos = [(1, -2), (0, -2), (-1, -2), (-2, -2),
           (1, -1), (0, -1), (-1, -1), (-2, -1)]
    pos += [(1, 0), (0, 0), (-1, 0), (-2, 0), (1, 1), (0, 1), (-1, 1), (-2, 1)]

    iccd = list(range(1, 17))

    return dict(zip(pos, iccd))


def align_quad(fields, params={}, j=0, output_q=None):

    print('Nfields', len(fields))

    df_obs = params['obs']
    ccd_map = map_ccd()
    df = pd.DataFrame()
    for field in fields:
        idx = df_obs['field'] == field
        sel = pd.DataFrame(df_obs[idx])
        align = Align_Field(sel, ccd_map=ccd_map)
        dfb = align()
        ndata_orig = len(sel)
        ndata_end = len(dfb)
        if ndata_end != ndata_orig:
            print('pb here', ndata_orig, ndata_end, field)
            plot_quad(sel, dfb)
        df = pd.concat((df, dfb))

    if output_q is not None:
        return output_q.put({j: df})
    else:
        return 0


def plot_quad(sel, dfb):

    import matplotlib.pyplot as plt
    rcids = dfb['rcid'].unique()
    print(rcids, len(rcids))
    fig, ax = plt.subplots()
    ax.plot(sel['ra'], sel['dec'], 'b.')
    ax.plot(dfb['ra'], dfb['dec'], 'ko', mfc='None')
    ax.plot(dfb['ra_quad'], dfb['dec_quad'], 'r*')
    """
    ax.plot(dfb['ra']*np.cos(np.deg2rad(dfb['dec'])), dfb['dec'], 'ko')
    ax.plot(dfb['ra_quad']*np.cos(np.deg2rad(dfb['dec_quad'])),
            dfb['dec_quad'], 'r*')
    """
    plt.show()
