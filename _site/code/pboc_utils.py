"""
This script contains an array of functions for use in the physical biology of
the cell course at the Gwangju Institute of Science and Technology in Gwangju,
PRK. These functions were written by Griffin Chure and carry an MIT license.

2017/03/10
MJM: removed skimage requirement for CSHL PBoC
    (Canopy w/ Python3 lacks it)

2017/03/15
MJM: reinstated skimage, phase_segmentation, and extract_intensities, 
    with slight variation from GC's original.
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import skimage.io
import skimage.segmentation
import skimage.filters
import skimage.measure
import skimage.morphology
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

def phase_segmentation(image, thresh, area_bounds=[1,3], ip_dist=0.16):
    """Take a phase contrast image, segment by thresholding, and return 
    the mask."""
    # first rescale image to intensities from 0 to 1
    im_float = (image - image.min()) / (image.max() - image.min())
    # do background subtraction
    im_blur = skimage.filters.gaussian(im_float, sigma=50.0)
    im_sub = im_float - im_blur
    # apply the threshold
    im_thresh = im_sub < thresh
    # next do area screen. but 1st need to label objects
    im_lab = skimage.measure.label(im_thresh)
    props = skimage.measure.regionprops(im_lab)
    # initialize an empty image
    approved_objects = np.zeros_like(im_lab)
    # loop over object in labeled image and apply area screen
    for labeled_obj in props:
        # convert pixel area to physical area in um^2
        area = labeled_obj.area * ip_dist**2
        # apply area screen
        if (area > area_bounds[0]) & (area < area_bounds[1]):
            approved_objects += (im_lab == labeled_obj.label)
    # clear border and relabel
    im_border = skimage.segmentation.clear_border(approved_objects)
    final_seg = skimage.measure.label(im_border)
    return final_seg

def extract_intensities(seg, fluo_im):
    """Takes two images as args: 1st is a segmentation mask, 
    2nd is a fluorescence image. Returns a list of intensities of each
    object in fluorescence image."""
    # get the regionprops of the fluorescent image, conditioned on the 
    # segmentation mask
    props = skimage.measure.regionprops(seg, intensity_image=fluo_im)
    cell_ints = []
    # loop over labeled objects
    for labeled_obj in props:
        # append intensity of each cell to list
        cell_ints.append(labeled_obj.mean_intensity)
    return np.array(cell_ints)