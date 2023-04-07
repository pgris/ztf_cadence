from ztf_pipeutils.ztf_util import multiproc
from ztf_metrics.metrics import CadenceMetric, RedMagMetric
import pandas as pd


def processMetric_multiproc(metricName, df, nproc, nside, coadd_night, npixels=-1, pixelList='all'):
    """
    Method to process metrics

    Parameters
    --------------
    metricName: str
      name of the metric
    df: pandas
      data to process
    nproc: int
      number of procs for multiprocessing
    nside: int
      nside healpix parameter
    coadd_night: int
      to perform coaddition of data per night and per band
    npixels: int,opt
      number of pixels to process (randomly chosen; default: -1: all)
    pixelList: str, opt
      list of pixels to process (default: all)

    Returns
    --------------
    New data frame for the differents pixels with the calculation (value) of differents metrics.
    """

    if pixelList != 'all':
        # set of pixels
        healpixIDs = pixelList.split(',')
    else:
        # all the pixels

        healpixIDs = ','.join(df['healpixID'].to_list())
        healpixIDs = set(healpixIDs.split(","))
        healpixIDs = list(filter(lambda a: a != 'pNonep', healpixIDs))

        """
        ll = df['healpixID'].to_list()
        import itertools
        healpixIDs = list(itertools.chain.from_iterable(ll))
        # healpixIDs = list(set(df['healpixID'].to_list()))
        """
    # random pixels
    if npixels >= 1:
        import random
        healpixIDs = random.choice(healpixIDs, k=npixels)

    # adjust nproc (if necessary)
    import numpy as np
    if nproc > np.min([nproc, len(healpixIDs)]):
        nproc = len(healpixIDs)

    print('finally', nproc, healpixIDs)
    params = {}
    params['metricName'] = metricName
    params['data'] = df
    params['nside'] = nside
    params['coadd_night'] = coadd_night

    resdf = multiproc(healpixIDs, params, processMetric, nproc)

    return resdf


def processMetric(healpixIDs, params={}, j=0, output_q=None):
    """
    Method to process metrics

    Parameters
    --------------
   healpixIDs: list
      list of healpixIDs to process
    params: dict, opt
      dict of parameters (default: {})

    Returns
    --------------
    New data frame for the differents pixels with the calculation (value) of differents metrics.
    """

    metricName = params['metricName']
    data = params['data']
    nside = params['nside']
    coadd_night = params['coadd_night']

    # print('pparams', params)

    healpixIDs = set(healpixIDs)

    cl = eval('{}(nside={},coadd_night={})'.format(
        metricName, nside, coadd_night))

    resdf = pd.DataFrame()
    fracs = range(10, 100, 50)
    prfr = dict(zip(fracs, [1]*len(fracs)))
    print('number of pixels to process', len(healpixIDs))
    #healpixIDs = ['p144577p']
    for hpix in healpixIDs:
        # print('processing', hpix, type(hpix))
        dfb = data[data['healpixID'].str.contains(hpix, regex=False)]

        df_new = dfb.copy()
        if len(df_new) < 10 and metricName != 'RedMagMetric':
            continue
        hpix = int(hpix.split('p')[1])
        respix = cl.run(int(hpix), df_new)
        resdf = pd.concat([resdf, respix])
        frac_processed = 100.*len(resdf)/len(healpixIDs)
        for key, vv in prfr.items():
            if frac_processed >= key and vv:
                print('fraction processed', frac_processed)
                prfr[key] = 0

    if output_q is not None:
        return output_q.put({j: resdf})
    else:
        return resdf
