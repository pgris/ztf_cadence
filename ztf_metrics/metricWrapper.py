from ztf_pipeutils.ztf_util import multiproc
from ztf_metrics.metrics import CadenceMetric
import pandas as pd


def processMetric_multiproc(metricName, df, nproc):
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

    healpixIDs = set(healpixIDs)

    cl = eval('{}()'.format(metricName))

    resdf = pd.DataFrame()
    for hpix in healpixIDs:
        dfb = data[data['healpixID'].str.contains(hpix)]
        df_new = dfb.copy()
        respix = cl.run(hpix, df_new)

        resdf = pd.concat([resdf, respix])
    if output_q is not None:
        return output_q.put({j: resdf})
    else:
        return resdf
