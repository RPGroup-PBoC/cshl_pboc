"""
This script contains an array of functions for use in the physical biology of
the cell course at the Gwangju Institute of Science and Technology in Gwangju,
PRK. These functions were written by Griffin Chure and carry an MIT license.

2017/03/10
MJM: removed skimage requirement for CSHL PBoC
    (Canopy w/ Python3 lacks it)
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# import skimage.io
# import skimage.segmentation
# import skimage.filters
# import skimage.measure
import seaborn as sns


def bar3(data, xlabel='x', ylabel='y', zlabel='z', bin_step=1,
         x_vec='default', y_vec='default'):
    """
    Generates a three-dimensional bar plot of provided data
    on a unique figure axis.

    Parameters
    ----------
    data : 2d-array
        Numpy array containing the z information for plotting. This should be
        NxM where N is the x axis and M is the y axis.
    xlabel : str
        Label for the x axis of the plot.
    ylabel : str
        Label for the y axis of the plot.
    zlabel : stk
        Label for the z axis of the plot.
    bin_step : int, default 1
        Resolution of plotting in the y dimension. For example, if '1' is
        provided, every y step will be plotted. If '5' is provided, every
        5th row will be plotted.
    x_vec : 1d-array, default is aranged.
        x-vector for plotting. If 'default', a linearly aranged vector is used.
    y_vec : 1d-array, default is aranged.
        y-vector for plotting. If 'default', a linearly aranged vector is used.

    Returns
    -------
    fig : matplotlib figure plotting axis.
        Figure object for plotting.

    ax  : matplotlib plotting axis.
        Axis object for further manipulation.
    """

    # Determine the x and y lengths.
    x_length, y_length = np.shape(data)

    # Instantiate the figure.
    fig = plt.figure()

    # Add a three-dimensional axis.
    ax = fig.add_subplot(1, 1, 1, projection='3d')

    # Set the colorscheme
    try:
        colors = sns.color_palette('viridis', n_colors=y_length)
    except:
        colors = sns.color_palette('RdBu_r', n_colors=y_length)

    # set up the plotting vectors.
    if x_vec is 'default':
        x_vec = np.arange(0, x_length, 1)
    if y_vec is 'default':
        y_vec = np.arange(0, y_length, 1)

    # Iterate through each y point and make the bar plot.
    for i in range(0, y_length, bin_step):
        ax.bar(x_vec, data[:, i], y_vec[i], zdir='x', width=1,
               color=colors[i])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    return fig, ax
