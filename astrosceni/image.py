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
  def __init__(self, path=None):
    self.path = ""
    self.original_data = None
    self.original_header = None
    self.cutout = None

    self.original_wcs = None  # Store WCS object for conversion
    self.cutout_wcs = None  # Store WCS object for conversion

    if (path != None): self.load(path)

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
    print("IMAGE_DATA", self.getImageData(original).shape)
    
    # coords = np.argwhere(data)
    # row_min, col_min = coords.min(axis=0)
    # row_max, col_max = coords.max(axis=0)

    row_min = 0
    col_min = 0
    row_max = data.shape[0] - 1
    col_max = data.shape[1] - 1

    if x_start == None: x_start = col_min
    if x_end == None: x_end = col_max
    if y_start == None: y_start = row_min
    if y_end == None: y_end = row_max
    
    print("COL_MAX", col_max)
    print("COL_MIN", col_min)

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

    Pass arguements in degrees or as string
    """
    data = self.getImageData(original=original)

    # Image pixel boundaries
    ny, nx = data.shape  # Shape: rows (Y) x columns (X)

    # Convert pixel boundaries to RA/Dec
    # Bottom-left (0, 0) and top-right (nx-1, ny-1)
    ra_min, dec_min = self.getWCS(original=original).all_pix2world(0, 0, 0)  # Bottom-left corner
    ra_max, dec_max = self.getWCS(original=original).all_pix2world(nx - 1, ny - 1, 0)  # Top-right corner

    if ra_start == None: ra_start = ra_min
    if ra_end == None: ra_end = ra_max
    if dec_start == None: dec_start = dec_min
    if dec_end == None: dec_end = dec_max

    # if isinstance(ra_start, str): ra_start = SkyCoord(ra=ra_start, dec="0d", unit=(u.hourangle, u.deg)).ra.deg
    # if isinstance(ra_end, str): ra_end = SkyCoord(ra=ra_end, dec="0d", unit=(u.hourangle, u.deg)).ra.deg
    # if isinstance(dec_start, str): dec_start = SkyCoord(ra="00h", dec=dec_start, unit=(u.hourangle, u.deg)).dec.deg
    # if isinstance(dec_end, str): dec_end = SkyCoord(ra="00h", dec=dec_end, unit=(u.hourangle, u.deg)).dec.deg
    ra_start = Image._convertDegRA(ra_start)
    ra_end = Image._convertDegRA(ra_end)
    dec_start = Image._convertDegDec(dec_start)
    dec_end = Image._convertDegDec(dec_end)

    # Convert RA/Dec to pixel coordinates using WCS
    position = self.getWCS(original=original).all_world2pix((ra_end+ra_start)/2, (dec_end+dec_start)/2, 0)
    position = (position[0], position[1])  # x, y in pixel coordinates
    pt1 = self.getWCS(original=original).all_world2pix(ra_start, dec_start, 0)
    pt2 = self.getWCS(original=original).all_world2pix(ra_end, dec_end, 0)
    size = (abs(pt2[1] - pt1[1]), abs(pt2[0] - pt1[0]))  # Height and width of the cutout

    # Create the cutout using Cutout2D
    cutout = Cutout2D(data, position, size, wcs=self.getWCS(original=original))

    # Update current image data and WCS
    self.cutout = cutout
    self.cutout_wcs = cutout.wcs  # Update WCS to match the cutout

  # def setLabeledStars(self, stars):
  #   for star in stars:


  # def applyContour(self, contour):

  # Assumes both images have exact same zoom
  # Consider zooming in for the lesser-zoomed in picture to fix this in the future?
  @staticmethod
  def overlap(img1, img2):
    # Return 2 new images that are the overlap of the image

    ra_min1, dec_min1, ra_max1, dec_max1 = img1.getBounds()
    ra_min2, dec_min2, ra_max2, dec_max2 = img2.getBounds()
    
    # print(dec_min1, dec_max1, dec_min2, dec_max2)
    is_ra_overlap, ra_overlap = Image._isOverlappedRA(ra_min1, ra_max1, ra_min2, ra_max2)
    is_dec_overlap, dec_overlap = Image._isOverlappedDec(dec_min1, dec_max1, dec_min2, dec_max2)

    if is_ra_overlap and is_dec_overlap:
      # img1.cropCoords(), img2
      return ra_overlap
    else: return False
    
  
  @staticmethod
  def _isOverlappedDec(start1, end1, start2, end2):
    """
    Returns whether two RA intervals overlap, arguements are in decimal degrees
    """
    start1 = Image._convertDegDec(start1)
    end1 = Image._convertDegDec(end1)
    start2 = Image._convertDegDec(start2)
    end2 = Image._convertDegDec(end2)

    # overlap = ((start2 > end1 and start2 < start1) or (start2 < end1 and start2 > start1)) or ((end2 > end1 and end2 < start1) or (end2 < end1 and start2 > start1))

    # return overlap

    # Ensure intervals are properly ordered
    start1, end1 = min(start1, end1), max(start1, end1)
    start2, end2 = min(start2, end2), max(start2, end2)

    # Check for overlap
    if max(start1, start2) <= min(end1, end2):  # Overlap condition
        overlap = True
        # Calculate the intersection interval
        intersection_start = max(start1, start2)
        intersection_end = min(end1, end2)
        return overlap, (intersection_start, intersection_end)
    else:
        overlap = False
        return overlap, None

  # @staticmethod
  # def _isOverlappedRA(start1, end1, start2, end2):
  #   """
  #   Returns whether two RA intervals overlap, arguements are in decimal degrees
  #   """
  #   start1 = Image._convertDegRA(start1)
  #   end1 = Image._convertDegRA(end1)
  #   start2 = Image._convertDegRA(start2)
  #   end2 = Image._convertDegRA(end2)

  #   overlap = ((start2 > end1 and start2 < start1) or (start2 < end1 and start2 > start1)) or ((end2 > end1 and end2 < start1) or (end2 < end1 and start2 > start1))

  #   return overlap

  @staticmethod
  def _isOverlappedDec(start1, end1, start2, end2):
    """
    Returns whether two Dec intervals overlap and the intersection interval if they overlap,
    considering the possibility of a reversed Dec axis.

    Arguments:
        start1, end1, start2, end2: Dec intervals in decimal degrees.
        wcs: WCS object to determine Dec orientation.

    Returns:
        overlap (bool): True if the intervals overlap, False otherwise.
        intersection (tuple): A tuple (start, end) representing the intersection interval in decimal degrees.
                              If there is no overlap, returns None.
    """
    # Detect Dec direction from the WCS object
    dec_reversed = wcs.wcs.cd[1, 1] < 0 if wcs.wcs.has_cd() else wcs.wcs.pc[1, 1] < 0

    # Reverse the intervals if Dec is reversed
    if dec_reversed:
        start1, end1 = -end1, -start1
        start2, end2 = -end2, -start2

    # Create SkyCoord objects for the Dec intervals
    dec1 = SkyCoord(ra=0 * u.deg, dec=[start1, end1] * u.deg)
    dec2 = SkyCoord(ra=0 * u.deg, dec=[start2, end2] * u.deg)

    # Sort the intervals to ensure start <= end
    start1, end1 = sorted(dec1.dec.degree)
    start2, end2 = sorted(dec2.dec.degree)

    # Check for overlap
    overlap = max(start1, start2) <= min(end1, end2)

    if overlap:
        # Calculate intersection
        intersection_start = max(start1, start2)
        intersection_end = min(end1, end2)

        # Reverse the intersection back if Dec is reversed
        if dec_reversed:
            intersection_start, intersection_end = -intersection_end, -intersection_start

        return overlap, (intersection_start, intersection_end)
    else:
        return overlap, None


  @staticmethod
  def _convertDegRA(ra):
    if isinstance(ra, str): return SkyCoord(ra=ra, dec="0d", unit=(u.hourangle, u.deg)).ra.deg
    return ra

  @staticmethod
  def _convertDegDec(dec):
    if isinstance(dec, str): return SkyCoord(ra="00h", dec=dec, unit=(u.hourangle, u.deg)).dec.deg
    return dec

  def plot(self, original=False, showCropped=False, croppedBorder='white', showLabeledStars=True, labelCircleSize=10):
    """
    Plots the cropped image by default, pass `original=True` to plot the original
    """
    data = self.getImageData(original)

    fig = plt.figure()

    # wcs = WCS(self.getWCS())
    ax = WCSAxes(fig, [0, 0, 1, 1], wcs=self.getWCS(original=original))
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

    if showCropped and self.cutout:
      self.cutout.plot_on_original(color=croppedBorder)

    # if self.labeled_stars is not None and showLabeledStars:
    #   self.labeled_stars.apply(lambda star: ax.add_patch(plt.Circle((star['x_pixels'], star['y_pixels']), labelCircleSize, color='b', fill=False)), axis=1)

    plt.show()

  def getWCS(self, original=False):
    if original: return self.original_wcs
    elif self.cutout: return self.cutout_wcs
    else: return self.original_wcs
  
  def getImageData(self, original=False):
    if original: return self.original_data
    elif self.cutout: return self.cutout.data
    else: return self.original_data
  
  def getBounds(self):
    """
    Returns the coordinates at the bounds of the cropped image
    """
    data = self.getImageData()
    ny, nx = data.shape  # Shape: rows (Y) x columns (X)

    ra_min, dec_min = self.getWCS().all_pix2world(0, 0, 0)  # Bottom-left corner
    ra_max, dec_max = self.getWCS().all_pix2world(nx - 1, ny - 1, 0)  # Top-right corner
    return ra_min, dec_min, ra_max, dec_max

  def plotHist(self, nbins=200):
    # flatten means: we put our 2d array in a 1d array
    histogram = plt.hist(self.getImageData().flatten(), nbins)

    plt.xlabel('Pixel Content')
    plt.ylabel('Number of Pixels')
    plt.yscale('log')
    plt.show()

  def setLabeledStars(self, stars_filter):
    self.labeled_stars = stars_filter.getVisibleStars()

  def _printInfo(self, original=False):
    if original: print("Original", self.original_data.shape)
    elif self.cutout != None: print("Original", self.cutout.data.shape)
    else: print("Original", self.original_data.shape)

  # Idea for overlap: find ra interval bounds in terms of pixels and then find the intersection interval in terms of pixels then convert back to ra_dec