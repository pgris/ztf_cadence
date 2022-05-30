import numpy as np
import pandas as pd


def loadData(dataDir):
    """
    Function to load data

    Parameters
    --------------
    dataDir: str
      directory where files are located

    Returns
    ----------
    pandas df of data

    """
    import glob
    fis = glob.glob('{}/*.hdf5'.format(dataDir))

    print(fis)
    data = pd.DataFrame()
    for fi in fis:
        tt = pd.read_hdf(fi)
        data = pd.concat((data, tt))

    return data


def binnedData(data, bins, xvar, yvar):
    """
    Function to build binned data

    Parameters
    --------------
    data: pandas df
       data to bin
    bins: list(float)
      bins used for the binning
    xvar: str
      x-axis variable
    yvar: str
      y-axis variable

    Returns
    -----------
    pandas df with binned data

    """
    bin_centers = (bins[: -1] + bins[1:])/2
    group = data.groupby(pd.cut(data[xvar], bins))
    bin_values = group.apply(
        lambda x: pd.DataFrame({yvar: [x[yvar].mean()], '{}_std'.format(yvar): [x[yvar].std()]}))

    bin_values[xvar] = bin_centers

    return bin_values
