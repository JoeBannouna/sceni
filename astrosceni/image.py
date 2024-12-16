import numpy as np
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
from matplotlib.colors import LogNorm, SymLogNorm
import os


class Image:
  def __init__(self):
    self.path = ""
    self.original_image_data = None
    self.current_image_data = None
    self.image_header = None

  def load(self, path):
    """
    Loads the image data and headers from a fits file
    """
    self.path = path
    filename = os.path.abspath(self.path)
    hdu = fits.open(filename)
    
    self.original_image_data = hdu[0].data
    self.current_image_data = hdu[0].data
    self.image_header = hdu[0].header

  def _printData(self, original=False):
    """
    Prints the pixel values of the cropped image, for original image pass `original=True`
    """
    if original: print(self.current_image_data)
    else: print(self.original_image_data)

  def cropPixels(self, x_start=None, x_end=None, y_start=None, y_end=None, original=True):
    """
    Creates a crop of the original image, to crop the currently cropped image pass `original=False`
    """
    data = []
    if original: data = self.original_image_data
    else: data = self.current_image_data
    
    coords = np.argwhere(data)
    row_min, col_min = coords.min(axis=0)
    row_max, col_max = coords.max(axis=0)

    if x_start == None: x_start = col_min
    if x_end == None: x_end = col_max
    if y_start == None: y_start = row_min
    if y_end == None: y_end = row_max

    print("x_start", x_start, "x_end", x_end, "y_start", y_start, "y_end", y_end)
    
    cropped = data[y_start:y_end+1, x_start:x_end+1]
    self.current_image_data = cropped

  # def cropCoords(self, contour):


  # def setLabeledStars(self, contour):

  # def applyContour(self, contour):

  def plot(self, original=False):
    """
    Plots the cropped image by default, pass `original=True` to plot the original
    """
    data = []
    if original: data = self.original_image_data
    else: data = self.current_image_data
    
    fig = plt.figure()

    wcs = WCS(self.image_header)
    ax = WCSAxes(fig, [0, 0, 1, 1], wcs=wcs)
    fig.add_axes(ax)

    # Now with an other colormap and in logscale
    img = ax.imshow(data, cmap='afmhot', norm=LogNorm(vmin=1, vmax=1e5))

    plt.scatter(0, 0)
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    ax.grid(color='white', ls='solid')

    # Set the limits of the colorbar and a label
    cb = fig.colorbar(img)
    # img.set_clim(1.1*np.min(image_data), np.max(image_data))
    cb.set_label('Counts')

    plt.show()
