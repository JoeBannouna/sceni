import copy
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
from matplotlib.colors import LogNorm, SymLogNorm
from astropy.nddata import Cutout2D
import os, warnings

#Simplifies warning message
warnings.simplefilter('always', UserWarning)
warnings.formatwarning = lambda message, *args: f"{message}\n"

class Image:
  def __init__(self, path=None):
    self.path = ""
    self.labeled_starts = None
    self.original_data = None
    self.original_header = None
    self.cutout = None
    self.cutout_data = np.array([])

    self.original_wcs = None  # Store WCS object for conversion
    self.cutout_wcs = None  # Store WCS object for conversion
    
    if (path != None): self.load(path)

  #Determines how many saturated pixels are in image, compares as ratio to resolution and provides warning if necessary
  def checkSaturatedPixelCount(self, saturatedRatioLim):
    count = np.sum(self.original_data > (2**15)-2)
    print("Count of saturated pixels: ", count)
    resolution = len(self.original_data) * len(self.original_data[0])
    print("Total image pixel count: ", resolution)

    ratio = count/resolution

    if (ratio >= saturatedRatioLim):
      warnings.warn("WARNING: Saturated pixels account for %.2f%% of total pixel count of image" % (100*ratio))

  def setSaturatedPixelsToNan(self, data):
    data[data > (2**15)-2] = np.nan
    return data
    

  def load(self, path, saturatedRatioLim = 0.001):
    """
    Loads the image data and headers from a fits file
    """
    self.path = path
    filename = os.path.abspath(self.path)
    hdu = fits.open(filename)
    
    self.original_data = hdu[0].data
    self.checkSaturatedPixelCount(saturatedRatioLim)
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
    self.cutout_data = cutout.data
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
    self.cutout_data = cutout.data
    self.cutout_wcs = cutout.wcs  # Update WCS to match the cutout

  @staticmethod
  def _convertDegRA(ra):
    """Converts 'ra' string/numerical input to degrees"""
    if isinstance(ra, str): return SkyCoord(ra=ra, dec="0d", unit=(u.hourangle, u.deg)).ra.deg
    return ra

  @staticmethod
  def _convertDegDec(dec):
    """Converts 'dec' string/numerical input to degrees"""
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

    if showCropped and self.cutout_data.size != 0:
      self.cutout.plot_on_original(color=croppedBorder)

    if self.labeled_starts is not None and showLabeledStars:
      self.labeled_starts.apply(lambda star: ax.add_patch(plt.Circle((star['x_pixels'], star['y_pixels']), labelCircleSize, color='b', fill=False)), axis=1)

    plt.show()

  def zoomToPoint(self, x = None, y = None, ra = None, dec = None, zoom_size = 50, original = False, custom_vmin = None, custom_vmax = None, custom_cmap = 'grey', show_colorbar = True,
                   setSaturatedToNan = False, use_log_norm = False):

    """
    Zooms into image around point
    - (x, y) pixel coordinates or
    - (ra, dec) world coordinates

    - zoom_size = half-width of zoom-in window in pixels
    - original = whether to zoom in from original image or already cropped image

    """

    data = self.getImageData(original)
    wcs = self.getWCS(original)

    if (setSaturatedToNan == True):
      data = self.setSaturatedPixelsToNan(data)
    print(np.nanmax(data))

    if ra is not None and dec is not None:
      x, y = wcs.all_world2pix(ra, dec, 0)

    if x is None or y is None:
      raise ValueError("Pixel coordinates (x, y) or world coordinates (ra, dec) must be provided")
    
    x, y = int(x), int(y)

    # Defining zoom window size
    x_start, x_end = max(0, x - zoom_size), min(x + zoom_size, data.shape[1])
    y_start, y_end = max(0, y - zoom_size), min(y + zoom_size, data.shape[0])

    #Extract zoom region
    zoomed_data = data[y_start:y_end, x_start:x_end]
    fig, ax = plt.subplots(figsize=(6, 6))

    if (custom_vmin == None):
      vmin = np.nanmin(data)
    else:
      vmin = custom_vmin
    if (custom_vmax == None):
      vmax = np.nanmax(data)
    else:
      vmax = custom_vmax

    norm = None
    if (use_log_norm == True):
      norm = LogNorm(vmin=vmin, vmax=vmax)
      img = ax.imshow(zoomed_data, origin = 'lower', cmap = custom_cmap, norm=norm)
    else:
      img = ax.imshow(zoomed_data, origin = 'lower', cmap = custom_cmap, vmin = vmin, vmax = vmax)

    ax.scatter(x-x_start, y-y_start, color = 'red', marker = '+', label = "Zoom Center")
    ax.set_xlim(0, x_end - x_start)
    ax.set_ylim(0, y_end - y_start)
    ax.set_title(f"Zoomed-in View (Center: x = {x}, y = {y})")
    ax.set_xlabel("X Pixel")
    ax.set_ylabel("Y Pixel")
    ax.legend()

    if (show_colorbar == True):
      cbar = fig.colorbar(img, ax = ax, fraction = 0.05, pad = 0.04)
      cbar.set_label("Pixel Intensity")

    plt.show()
    

  def getWCS(self, original=False):
    if original: return self.original_wcs
    elif self.cutout_data.size != 0: return self.cutout_wcs
    else: return self.original_wcs
    
  def setWCS(self, wcs, original=False):
    if original: self.original_wcs = wcs
    elif self.cutout_data.size != 0: self.cutout_wcs = wcs
    else: self.original_wcs = wcs
  
  def getImageData(self, original=False):
    if original: return self.original_data
    elif self.cutout_data.size != 0: return self.cutout_data
    else: return self.original_data
  
  def setImageData(self, data, original=False):
    """
    When calling this function its important to remember that WCS is invalid unless setWCS is called along with this function to ensure the data is coherent
    """
    if original: self.original_data = data
    elif self.cutout_data.size != 0: self.cutout_data = data
    else: self.original_data = data
  
  def getBounds(self):
    """
    Returns the coordinates at the bounds of the cropped image
    """
    data = self.getImageData()
    ny, nx = data.shape  # Shape: rows (Y) x columns (X)

    ra_min, dec_min = self.getWCS().all_pix2world(0, 0, 0)  # Bottom-left corner
    ra_max, dec_max = self.getWCS().all_pix2world(nx - 1, ny - 1, 0)  # Top-right corner
    return ra_min, dec_min, ra_max, dec_max

  @staticmethod
  def subtract(NB_image, BB_image, mu=1):
    result = copy.deepcopy(NB_image)
    
    nbdata = np.array([])
    if NB_image.cutout != None: nbdata = NB_image.cutout_data
    else: nbdata = NB_image.original_data

    bbdata = np.array([])
    if BB_image.cutout != None: bbdata = BB_image.cutout_data
    else: bbdata = BB_image.original_data

    result.cutout_data = nbdata - mu*bbdata
    result.cutout_wcs = NB_image.getWCS()

    return result

  def setLabeledStars(self, stars_filter):
    self.labeled_starts = stars_filter.getVisibleStars()

  # def applyContour(self, contour):

  def plotHist(self, nbins=200):
    # flatten means: we put our 2d array in a 1d array
    histogram = plt.hist(self.getImageData().flatten(), nbins)

    plt.xlabel('Pixel Content')
    plt.ylabel('Number of Pixels')
    plt.yscale('log')
    plt.show()