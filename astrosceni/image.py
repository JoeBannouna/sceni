import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
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

  def cropCoords(self, ra_start=None, ra_end=None, dec_start=None, dec_end=None, original=True):
    """
    Crops the image based on RA and Dec coordinates.
    If None is passed, defaults to the image's RA/Dec boundaries.

    ra_start, ra_end, dec_start, dec_end: RA/Dec coordinates in degrees.
    """
    data = self.getImageData(original)

    # Image pixel boundaries
    ny, nx = data.shape  # Shape: rows (Y) x columns (X)

    # Convert pixel boundaries to RA/Dec
    # Bottom-left (0, 0) and top-right (nx-1, ny-1)
    ra_min, dec_min = self.original_wcs.all_pix2world(0, 0, 0)  # Bottom-left corner
    ra_max, dec_max = self.original_wcs.all_pix2world(nx - 1, ny - 1, 0)  # Top-right corner

    c = SkyCoord('00h42m30s', '+41d12m00s', frame='icrs')

    if ra_start == None: ra_start = ra_min
    if ra_end == None: ra_end = ra_max
    if dec_start == None: dec_start = dec_min
    if dec_end == None: dec_end = dec_max

    if isinstance(ra_start, str): ra_start = SkyCoord(ra=ra_start, dec="0d", unit=(u.hourangle, u.deg)).ra.deg
    if isinstance(ra_end, str): ra_end = SkyCoord(ra=ra_end, dec="0d", unit=(u.hourangle, u.deg)).ra.deg
    if isinstance(dec_start, str): dec_start = SkyCoord(ra="00h", dec=dec_start, unit=(u.hourangle, u.deg)).dec.deg
    if isinstance(dec_end, str): dec_end = SkyCoord(ra="00h", dec=dec_end, unit=(u.hourangle, u.deg)).dec.deg

    print("ra_start", ra_start, "ra_end", ra_end, "dec_start", dec_start, "dec_end", dec_end)

    # Convert RA/Dec to pixel coordinates using WCS
    position = self.original_wcs.all_world2pix((ra_end+ra_start)/2, (dec_end+dec_start)/2, 0)
    position = (position[0], position[1])  # x, y in pixel coordinates
    pt1 = self.original_wcs.all_world2pix(ra_start, dec_start, 0)
    pt2 = self.original_wcs.all_world2pix(ra_end, dec_end, 0)
    size = (abs(pt2[1] - pt1[1]), abs(pt2[0] - pt1[0]))  # Height and width of the cutout

    print(position)
    print(size)

    # Create the cutout using Cutout2D
    cutout = Cutout2D(data, position, size, wcs=self.original_wcs)

    # Update current image data and WCS
    self.cutout = cutout
    self.cutout_wcs = cutout.wcs  # Update WCS to match the cutout

  # def setLabeledStars(self, contour):

  # def applyContour(self, contour):

  def plot(self, original=False):
    """
    Plots the cropped image by default, pass `original=True` to plot the original
    """
    data = self.getImageData(original)
    
    fig = plt.figure()

    # wcs = WCS(self.getWCS())
    ax = WCSAxes(fig, [0, 0, 1, 1], wcs=self.getWCS())
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
  
  def getImageData(self, original=False):
    if original: return self.original_data
    elif self.cutout: return self.cutout.data
    else: return self.original_data