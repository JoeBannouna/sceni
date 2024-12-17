import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
from matplotlib.colors import LogNorm, SymLogNorm
from astropy.nddata import Cutout2D

import os

class Image:
  def __init__(self):
    self.path = ""
    self.original_data = None
    self.original_header = None
    self.cutout = None

    self.original_wcs = None  # Store WCS object for conversion
    self.cutout_wcs = None  # Store WCS object for conversion

  def load(self, path):
    """
    Loads the image data and headers from a fits file
    """
    self.path = path
    filename = os.path.abspath(self.path)
    hdu = fits.open(filename)
    
    self.original_data = hdu[0].data
    self.original_header = hdu[0].header

    # Initialize WCS object from the FITS header
    self.original_wcs = WCS(self.original_header)

  def _printData(self, original=False):
    """
    Prints the pixel values of the cropped image, for original image pass `original=True`
    """
    print(self.getImageData(original))

  def cropPixels(self, x_start=None, x_end=None, y_start=None, y_end=None, original=True):
    """
    Creates a cutout of the original image, to crop the currently cropped image pass `original=False`
    """
    data = self.getImageData(original)
    
    coords = np.argwhere(data)
    row_min, col_min = coords.min(axis=0)
    row_max, col_max = coords.max(axis=0)

    if x_start == None: x_start = col_min
    if x_end == None: x_end = col_max
    if y_start == None: y_start = row_min
    if y_end == None: y_end = row_max
    
    if x_end < 0: x_end = col_max + x_end
    if y_end < 0: y_end = row_max + y_end

    position = ((x_start+x_end)/2, (y_start+y_end)/2)
    size = (y_end-y_start, x_end-x_start)
    cutout = Cutout2D(data, position, size, wcs=self.original_wcs)
    
    self.cutout = cutout
    self.cutout_wcs = cutout.wcs

  # def cropCoords(self, ra_start=None, ra_end=None, dec_start=None, dec_end=None, original=True):
  #     """
  #     Crops the image based on RA and Dec coordinates.
  #     ra_start, ra_end, dec_start, dec_end: RA/Dec coordinates in degrees.
  #     """
  #     data = []
  #     if original:
  #         data = self.original_data
  #     else:
  #         data = self.current_data



  #     # Convert RA/Dec to pixel coordinates
  #     # The WCS object can convert world coordinates (RA/Dec) to pixel coordinates
  #     x_start, y_start = self.wcs.all_world2pix(ra_start, dec_start, 0)
  #     x_end, y_end = self.wcs.all_world2pix(ra_end, dec_end, 0)


  #     # x_start = int(max(0, x_start)) if ra_start == None else None

  #     # Ensure the values are within the bounds of the image
  #     x_start, x_end = int(max(0, x_start)), int(min(data.shape[1] - 1, x_end))
  #     y_start, y_end = int(max(0, y_start)), int(min(data.shape[0] - 1, y_end))

  #     # Crop the data based on pixel indices calculated from RA/Dec
  #     cropped = data[y_start:y_end+1, x_start:x_end+1]
  #     self.current_data = cropped

  # def cropCoords(self, ra_start=None, ra_end=None, dec_start=None, dec_end=None, original=True):
  #   """
  #   Crops the image based on RA and Dec coordinates.
  #   If None is passed, defaults to the image's RA/Dec boundaries.

  #   ra_start, ra_end, dec_start, dec_end: RA/Dec coordinates in degrees.
  #   """
  #   data = self.original_data if original else self.current_data

  #   # Image pixel boundaries
  #   ny, nx = data.shape  # Shape: rows (Y) x columns (X)

  #   # Convert pixel boundaries to RA/Dec
  #   # Bottom-left (0, 0) and top-right (nx-1, ny-1)
  #   ra_bl, dec_bl = self.wcs.all_pix2world(0, 0, 0)  # Bottom-left corner
  #   ra_tr, dec_tr = self.wcs.all_pix2world(nx - 1, ny - 1, 0)  # Top-right corner

  #   # Assign defaults if None
  #   if ra_start is None:
  #       ra_start = min(ra_bl, ra_tr)
  #   if ra_end is None:
  #       ra_end = max(ra_bl, ra_tr)
  #   if dec_start is None:
  #       dec_start = min(dec_bl, dec_tr)
  #   if dec_end is None:
  #       dec_end = max(dec_bl, dec_tr)

  #   # Convert RA/Dec to pixel coordinates
  #   x_start, y_start = self.wcs.all_world2pix(ra_start, dec_start, 0)
  #   x_end, y_end = self.wcs.all_world2pix(ra_end, dec_end, 0)

  #   # Ensure the pixel coordinates are within bounds
  #   x_start, x_end = int(max(0, x_start)), int(min(nx - 1, x_end))
  #   y_start, y_end = int(max(0, y_start)), int(min(ny - 1, y_end))

  #   # Crop the data based on pixel indices
  #   cropped = data[y_start:y_end+1, x_start:x_end+1]
  #   self.current_data = cropped
  #   self.current_header.update(self.wcs[y_start:y_end+1, x_start:x_end+1].to_header())


  # def cropCoords(self, ra_start=None, ra_end=None, dec_start=None, dec_end=None, original=True):
  #   """
  #   Crops the image based on RA and Dec coordinates.
  #   If None is passed, defaults to the image's RA/Dec boundaries.

  #   ra_start, ra_end, dec_start, dec_end: RA/Dec coordinates in degrees.
  #   """
  #   data = self.getImageData(original)

  #   # Image pixel boundaries
  #   ny, nx = data.shape  # Shape: rows (Y) x columns (X)

  #   # Convert pixel boundaries to RA/Dec
  #   # Bottom-left (0, 0) and top-right (nx-1, ny-1)
  #   ra_bl, dec_bl = self.wcs.all_pix2world(0, 0, 0)  # Bottom-left corner
  #   ra_tr, dec_tr = self.wcs.all_pix2world(nx - 1, ny - 1, 0)  # Top-right corner

  #   # Assign defaults if None
  #   if ra_start is None:
  #       ra_start = min(ra_bl, ra_tr)
  #   if ra_end is None:
  #       ra_end = max(ra_bl, ra_tr)
  #   if dec_start is None:
  #       dec_start = min(dec_bl, dec_tr)
  #   if dec_end is None:
  #       dec_end = max(dec_bl, dec_tr)

  #   # Convert RA/Dec to pixel coordinates
  #   x_start, y_start = self.wcs.all_world2pix(ra_start, dec_start, 0)
  #   x_end, y_end = self.wcs.all_world2pix(ra_end, dec_end, 0)

  #   # Ensure the pixel coordinates are within bounds
  #   x_start, x_end = int(max(0, x_start)), int(min(nx - 1, x_end))
  #   y_start, y_end = int(max(0, y_start)), int(min(ny - 1, y_end))

  #   # Crop the data based on pixel indices
  #   cropped = data[y_start:y_end+1, x_start:x_end+1]
  #   self.current_data = cropped
  #   self.current_header.update(self.wcs[y_start:y_end+1, x_start:x_end+1].to_header())




  # def cropCoords(self, contour):


  # def setLabeledStars(self, contour):

  # def applyContour(self, contour):

  def plot(self, original=False):
    """
    Plots the cropped image by default, pass `original=True` to plot the original
    """
    data = self.getImageData(original)
    
    fig = plt.figure()

    # wcs = WCS(self.getCurrentWCS())
    ax = WCSAxes(fig, [0, 0, 1, 1], wcs=self.getCurrentWCS())
    fig.add_axes(ax)

    # Now with an other colormap and in logscale
    img = ax.imshow(data, cmap='afmhot', norm=LogNorm())

    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    ax.grid(color='white', ls='solid')

    # Set the limits of the colorbar and a label
    cb = fig.colorbar(img)
    # img.set_clim(1.1*np.min(image_data), np.max(image_data))
    cb.set_label('Counts')

    plt.show()

  def getWCS(self, orignial=False):
    if self.cutout: return self.cutout_wcs
    return self.original_wcs
  
  def getImageData(self, original):
    if original: return self.original_data
    elif self.cutout: return self.cutout.data
    else: return self.original_data