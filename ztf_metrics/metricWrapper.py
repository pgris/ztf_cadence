from ztf_pipeutils.ztf_util import multiproc
from ztf_metrics.metrics import CadenceMetric
import pandas as pd


def processMetric_multiproc(metricName, df, nproc, nside, coadd_night):
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

    Returns
    --------------
    New data frame for the differents pixels with the calculation (value) of differents metrics.
    """

    healpixIDs = ','.join(df['healpixID'].to_list())
    healpixIDs = healpixIDs.split(",")
    if 'None' in healpixIDs:
        healpixIDs.remove('None')

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

    healpixIDs = set(healpixIDs)

    cl = eval('{}(nside={},coadd_night={})'.format(
        metricName, nside, coadd_night))

    resdf = pd.DataFrame()
    for hpix in healpixIDs:
        dfb = data[data['healpixID'].str.contains(hpix)]
        df_new = dfb.copy()
        respix = cl.run(int(hpix), df_new)

        resdf = pd.concat([resdf, respix])

    if output_q is not None:
        return output_q.put({j: resdf})
    else:
        return resdf
