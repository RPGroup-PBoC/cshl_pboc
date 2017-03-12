"""
This script contains an array of functions for use in the physical biology of
the cell course at the Gwangju Institute of Science and Technology in Gwangju,
PRK. These functions were written by Griffin Chure and carry an MIT license.

2017/03/10
MJM: removed skimage requirement for CSHL PBoC
    (Canopy w/ Python3.5 lacks it)
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# import skimage.io
# import skimage.segmentation
# import skimage.filters
# import skimage.measure
import seaborn as sns

# # We'll start with the segmentation function.
# def phase_segmentation(image, threshold, area_bounds=[0.5, 6.0],
#                        ip_dist=0.160):
#     """
#     Segement a phase image and return the mask.
#
#     Parameters
#     ----------
#     image : 2d-array
#         The phase image to be segmented. This image will be converted to a
#         float type.
#     threshold : float
#         Threshold value for the segmentation. This function will select objects
#         below this threshold value.
#     area_bounds : list, default=[0.5, 6.0]
#         Area bounds for identified objects. This should be a list of two entries.
#     ip_dist : int or float, default = 0.160
#         Interpixel distance for the camera. This should be in units of microns
#         per pixel.
#
#     Returns
#     -------
#     final_seg : 2d-array
#         Final, labeled segmentation mask.
#     """
#
#     # First is to convert the image to a float.
#     im_float = (image - image.min()) / (image.max() - image.min())
#
#     # Do a background subtraction.
#     im_blur = skimage.filters.gaussian(im_float, sigma=50.0)
#     im_sub = im_float - im_blur
#
#     # Apply the threshold.
#     im_thresh = im_sub < threshold  # Note that we are using the provided arg
#
#     # Label the image and apply the area bounds.
#     im_lab = skimage.measure.label(im_thresh)
#     props = skimage.measure.regionprops(im_lab)
#     approved_objects = np.zeros_like(im_lab)
#     for prop in props:
#         area = prop.area * ip_dist**2
#         if (area > area_bounds[0]) & (area < area_bounds[1]):
#             approved_objects += im_lab == prop.label
#
#     # Clear the border and relabel.
#     im_border = skimage.segmentation.clear_border(approved_objects > 0)
#     final_seg = skimage.measure.label(im_border)
#
#     # Return the final segmentation mask
#     return final_seg
#
#
# # Now let's try writing one for to extract the mean intensities.
# def extract_intensity(seg, fluo_im):
#     """
#     Extracts the mean intensity of objects in a segmented image.
#
#     Parameters
#     ----------
#     seg : 2d-array, int
#         Segmentation mask with labeled objects.
#     fluo_im : 2d-array, int
#         Fluorescence image to extract intensities from.
#
#     Returns
#     -------
#     cell_ints : 1d-array
#         Vector of mean cell intensities. This has a length the same as the
#         number of objects in the provided segmentation mask.
#     """
#
#     # Get the region props of the fluorescence image using the segmentation
#     # mask.
#     props = skimage.measure.regionprops(seg, intensity_image=fluo_im)
#     cell_ints = []
#     for prop in props:
#         cell_ints.append(prop.mean_intensity)
#
#     # Convert the cell_ints to an array and return.
#     return np.array(cell_ints)


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
