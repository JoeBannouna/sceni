import copy
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.visualization import simple_norm
from astropy.visualization.wcsaxes import WCSAxes
from matplotlib.colors import LogNorm, SymLogNorm
from astropy.nddata import Cutout2D
from scipy.ndimage import gaussian_filter
import os, warnings

#Simplifies warning message
warnings.simplefilter('always', UserWarning)
warnings.formatwarning = lambda message, *args: f"{message}\n"
# Above 2 lines of warning generation code comes from the astropy library.


class Image:
  def __init__(self, path=None):
    self.path = ""
    self.labeled_starts = None
    self.original_data = None
    self.original_header = None
    self.cutout = None
    self.cutout_data = np.array([])
    self.saturated_pixels_removed = False

    self.original_wcs = None  # Store WCS object for conversion
    self.cutout_wcs = None  # Store WCS object for conversion
    
    if (path != None): self.load(path)


  #Determines how many saturated pixels are in image, compares as ratio to resolution and provides warning if necessary
  def checkSaturatedPixelCount(self, saturated_ratio_lim, saturation_threshold=None):
    """
    Determines the number of saturated pixels in an image and compares it as a ratio to the image resolution. A warning will be produced if necessary.

    Parameters
    ----------
    saturated_ratio_lim: float
      If the ratio (<saturated pixel counts>/<resolution>) is greater than or equal to saturated ratio limit, a warning will be produced.
      Resolution is the total number of image pixels.

    saturation_threshold: float
      Optional.
      Default: The maximum pixel value.
      The maximum pixel value.
      Each image has a maximum pixel value (max brightness) for which any greater values (brightness) will be set to this max pixel value.
      
    Returns
    -------
    count: integer
      The number of saturated pixels.
    """
    if self.saturated_pixels_removed:
      return 0
    if saturation_threshold == None: saturation_threshold = np.max(self.getImageData())
    count = np.sum(self.getImageData() >= saturation_threshold)
    resolution = len(self.getImageData()) * len(self.getImageData()[0]) # total no. pixels in the image
    ratio = count/resolution

    if (ratio >= saturated_ratio_lim):
      warnings.warn("WARNING: Saturated pixels account for %.2f%% of total pixel count of image" % (100*ratio))
    return count


  def setSaturatedPixelsToNan(self, saturation_threshold=None):
    """
    Changes the pixel values of saturated pixels to NaN.

    Parameter
    ---------
    saturation_threshold: float
      Optional
      Default: Maximum pixel value.
      Sets the value for pixels to be considered saturated.
    """

    if self.saturated_pixels_removed: return
    if saturation_threshold == None: saturation_threshold = np.max(self.getImageData())
    data = self.getImageData()
    data[data >= saturation_threshold] = np.nan
    self.setImageData(data)
    self.saturated_pixels_removed = True


  def load(self, path, saturated_ratio_lim = 0.001):
    """
    Loads a given image data and headers from a fits file.

    Parameters
    ----------
    path: string
      Path to the image file from the directory "SCENI".

    saturated_ratio_lim: float
      Optional
      Default: 0.001
      If the ratio (<saturated pixel counts>/<resolution>) is greater than or equal to saturated ratio limit, a warning will be produced.
      Resolution is the total number of image pixels.
    """
    self.path = path
    filename = os.path.abspath(self.path)
    hdu = fits.open(filename)
    
    self.original_data = hdu[0].data
    self.original_header = hdu[0].header
    self.checkSaturatedPixelCount(saturated_ratio_lim)

    # Initialize WCS object from the FITS header
    self.original_wcs = WCS(self.original_header)


  def cropPixels(self, x_start=None, x_end=None, y_start=None, y_end=None, original=True):
    """
    Creates a cutout of the original image, to crop the currently cropped image pass `original=False`
    """
    data = self.getImageData(original)
    ny, nx = data.shape

    if x_start == None: x_start = 0
    if x_end == None: x_end = nx #  Set x_end to last x pixel value if no value is passed
    if y_start == None: y_start = 0 
    if y_end == None: y_end = ny # Set y_end to last y pixel value if no value is passed
    
    # If negative values for x_end and y_end are passed, go backwards from the ends of the image
    if x_end < 0: x_end = nx + x_end 
    if y_end < 0: y_end = ny + y_end
    

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
    position = self.getWCS(original=original).all_world2pix((ra_end+ra_start)/2, (dec_end+dec_start)/2, 0) # The cutout array's center with respect to the data array
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
    """
    Converts string/numerical input for ra (right ascension) to degrees.

    Parameter
    ----------
    ra: string, float
      The right ascension coordinate.
    """
    if isinstance(ra, str): return SkyCoord(ra=ra, dec="0d", unit=(u.hourangle, u.deg)).ra.deg
    return ra


  @staticmethod
  def _convertDegDec(dec):
    """
    Converts string/numerical input for dec (declination) to degrees.

    Parameter
    ---------
    dec: string, float
      The declination coordinate.
    """
    if isinstance(dec, str): return SkyCoord(ra="00h", dec=dec, unit=(u.hourangle, u.deg)).dec.deg
    return dec


  def plot(self, original=False, showCropped=False, croppedBorder='white', showLabeledStars=False, labelCircleSize=10, labelCircleColor='b', cmap='afmhot', mode='linear', export_path = None):
    """
    Plots the cropped image by default, pass `original=True` to plot the original

    Parameters
    ----------
    original: boolean
      Optional
      Default: False (i.e. plot the cropped image)
      Controls whether or not to plot the original image.
      Pass 'original=True' to plot the original image.

    showCropped: boolean
      Optional
      Default: False
      Controls whether or not to show the cropped region.
      If showCropped = True, the original image will be plotted and the cropped region will be overlayed and indicated by a box.
    
    croppedBorder: string
      Optional
      Default: 'white'
      Sets the color of the box indicating the cropped region.

    showLabeledStars: boolean
      Optional
      Default: False
      Controls whether or not to show labelled stars.
      Labelled stars will be surrounded by a colored circle.

    labelCircleSize: integer
      Optional
      Default: 10
      Sets the size of the circle surrounding labelled stars.

    labelCircleColor: string
      Optional
      Default: 'b' (blue)
      Sets the color of the circle surrounding labelled stars.

    export_path: string
      Optional
      Default: None
      Sets the string to be where the plot is exported. Doesn't export if no string given
    """
    if showCropped: original = True
    data = self.getImageData(original)
    
    fig = plt.figure()

    # wcs = WCS(self.getWCS())
    ax = WCSAxes(fig, [0, 0, 1, 1], wcs=self.getWCS(original=original))
    fig.add_axes(ax)

    # Now with an other colormap and in logscale
    img = ax.imshow(data, cmap=cmap, norm=None if mode != 'log' else LogNorm(), 
                    vmin=None if mode != 'linear' else np.nanpercentile(self.getImageData(), 5), 
                    vmax=None if mode != 'linear' else np.nanpercentile(self.getImageData(), 99)
                    )

    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    ax.grid(color='white', ls='solid')

    # Set the limits of the colorbar and a label
    cb = fig.colorbar(img)
    # img.set_clim(1.1*np.min(image_data), np.max(image_data))
    cb.set_label('Counts')

    if showCropped and self.cutout_data.size != 0:
      print(self.cutout_data.size)
      self.cutout.plot_on_original(color=croppedBorder)

    if self.labeled_starts is not None and showLabeledStars:
      self.labeled_starts.apply(lambda star: ax.add_patch(plt.Circle(SkyCoord(star['RA'], star['DEC'], unit='deg').to_pixel(self.original_wcs if original else self.getWCS()), labelCircleSize, color=labelCircleColor, fill=False)), axis=1)

    if export_path:
      fig.savefig(export_path, bbox_inches = 'tight', pad_inches = 0.2)

    plt.show()


  def zoomToPoint(self, x = None, y = None, ra = None, dec = None, zoom_size = 50, original = False, custom_vmin = None, custom_vmax = None, custom_cmap = 'grey', show_colorbar = True,
                   setSaturatedToNan = False, use_log_norm = False):

    """
    Zooms into image around point
    - (x, y) pixel coordinates or
    - (ra, dec) world coordinates

    Parameters
    ----------
    x: integer
      Pixel coordinate along the x-axis.

    y: integer
      Pixel coordinate along the y-axis.

    ra: float
      Right ascension coordinate.

    dec: float
      Declination coordinate.

    zoom_size: integer
      Optional
      Default: 50
      Sets the magnification. It is the half-width of the zoom-in window in pixels.

    original: boolean
      Optional
      Default: False (i.e. zooms from the cropped image)
      Controls whether to zoom from the original image (if True) or from the cropped image (if False).

    custom_vmin: integer
      Optional
      Default: Minimum pixel value
      Sets the value for the minimum value of the colorbar when plotting the image.

    custom_vmax: integer
      Optional
      Default: Maximum pixel value
      Sets the value for the maximum value of the colorbar when plotting the image.

    custom_cmap: string
      Optional
      Default: grey
      Sets the colormap when plotting the image.

    show_colorbar: boolean
      Optional
      Default: True
      Controls whether or not to show the colorbar.

    setSaturatedToNan: boolean
      Optional 
      Default: False
      Controls whether or not to set the pixel values of saturated pixels to NaN.
      Each image has a maximum pixel value (max brightness) for which any greater values (brightness) will be set to this max pixel value. Pixels with such a value are saturated.
    
    use_log_norm: boolean
      Optional
      Default: False
      Controls whether or not to plot the image using LogNorm.
    """

    data = self.getImageData(original)
    wcs = self.getWCS(original)

    if (setSaturatedToNan == True):
      data = self.setSaturatedPixelsToNan(data)

    if ra is not None and dec is not None:
      x, y = wcs.all_world2pix(ra, dec, 0) # convert ra and dec into pixel coords, 0 is the original parameter which specifies the reference point

    if x is None or y is None:
      raise ValueError("Pixel coordinates (x, y) or world coordinates (ra, dec) must be provided")
    
    x, y = int(x), int(y)

    # Defining zoom window size
    x_start, x_end = max(0, x - zoom_size), min(x + zoom_size, data.shape[1]) # Restrict x_start >= 0; x - zoom_size is the half width of the zoomed image along x-axis
    y_start, y_end = max(0, y - zoom_size), min(y + zoom_size, data.shape[0]) # Restrict y_start >= 0; y - zoom_size is the half width of the zoom image along y-axis

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
    """
    Gets the header information from FITS files.;

    Parameters:
    -----------
    original: boolean
      Optional
      Default: False
      Controls whether or not to get the header information from the original FITS file. Otherwise, the cropped FITS file will be used.
    """
    
    if original: return self.original_wcs
    elif self.cutout_data.size != 0: return self.cutout_wcs
    else: return self.original_wcs
    

  def setWCS(self, wcs, original=False):
    """
    Sets header information to a FITS file.

    Parameters:
    -----------
    wcs: 
      Header information.
    
    original: boolean
      Optional
      Default: False
      Controls whether to add header information to the original or cropped image.
    """

    if original: self.original_wcs = wcs
    elif self.cutout_data.size != 0: self.cutout_wcs = wcs
    else: self.original_wcs = wcs
  

  def getImageData(self, original=False):
    """
    Gets the pixel data from an image.

    Parameter
    ---------
    original: boolean
      Optional
      Default: False
      Controls whether to use the original (if True) or cropped image (if False).
    """

    if original: return self.original_data
    elif self.cutout_data.size != 0: return self.cutout_data
    else: return self.original_data
  

  def setImageData(self, data, original=False):
    """
    Sets the pixel data for an image.
    When calling this function its important to remember that WCS is invalid unless setWCS is called along with this function to ensure the data is coherent.

    Parameters
    ----------
    data: numpy array
      The pixel data.
    
    original: boolean
      Optional
      Default: False
      Controls whether to set the pixel data to the original (if True) or cropped image (if False).
    """
    if original: self.original_data = data
    elif self.cutout_data.size != 0: self.cutout_data = data
    else: self.original_data = data
  

  def getBounds(self):
    """
    Returns the ra and dec coordinates at the bounds of the cropped image.

    Returns
    -------
    Ra and dec coordinates at the bounds of the cropped image.
    """
    data = self.getImageData()
    ny, nx = data.shape  # Shape: rows (Y) x columns (X)

    ra_bl, dec_bl = self.getWCS().all_pix2world(0, 0, 0)  # Bottom-left corner
    ra_br, dec_br = self.getWCS().all_pix2world(nx - 1, 0, 0)  # Bottom-right corner
    ra_tl, dec_tl = self.getWCS().all_pix2world(0, ny - 1, 0)  # Top-left corner
    ra_tr, dec_tr = self.getWCS().all_pix2world(nx - 1, ny - 1, 0)  # Top-right corner
    
    return SkyCoord(ra=ra_bl, dec=dec_bl, unit='deg'), SkyCoord(ra=ra_br, dec=dec_br, unit='deg'), SkyCoord(ra=ra_tl, dec=dec_tl, unit='deg'), SkyCoord(ra=ra_tr, dec=dec_tr, unit='deg')


  @staticmethod
  def subtract(NB_image, BB_image, mu=1):
    """
    Subtracts the pixel values between the NB and BB image according to the equation:
      NB - mu*BB
    
    Parameters:
    NB_Image: numpy array
      Numpy array of the pixel values for the NB (narrowband) image.

    BB_Image: numpy array
      Numpy array of the pixel values for the BB (broadband image).

    mu: float
      Optional
      Default: 1
      Scaling parameter for the subtraction.

    Returns
    -------
    result: numpy array
      The result of NB - mu*BB.
    """

    if NB_image.getImageData().shape != BB_image.getImageData().shape: raise ValueError('Images must be the same size')

    result = copy.deepcopy(NB_image)
    result.setImageData(NB_image.getImageData() - mu*BB_image.getImageData())
    result.saturated_pixels_removed = True 
    result.labeled_starts = None
    return result


  def setLabeledStars(self, stars_filter):
    """
    Gets a list of visible stars within a specified region.

    Parameter
    ---------
    stars_filter: pandas dataframe
      Visible stars within the region.
    """

    self.labeled_starts = stars_filter.getVisibleStars()


  def plotContour(self, sigma = None, levels = 5, cmap = 'viridis', norm_type = 'linear', base_cmap = 'gray', overlay = False, alpha = 0.5):
    """
    Plots a contour over the image which outlines the boundaries of all features.

    Parameters
    ----------
    sigma: integer
      Determines how blurry the image is.
      The higher the sigma, the cleaner the contours.

    levels: integer
      Optional
      Default: 5
      Number of contours to plot.

    cmap: string
      Optional
      Default: 'viridis'
      Sets the colormap for the contour plot.

    norm_type: string
      Optional
      Default: 'linear'
      Type of normalisation to use for the colormap of the contour.


    base_cmap: string
      Optional
      Default: 'gray'
      Sets the colormap of the base image.

    overlay: boolean
      Optional
      Default: False
      Controls whether or not to overlay the contour over the plotted image. 
      If overlay is False, only the contour is plotted, the image will not be plotted.

    alpha: float
      Optional
      Default: 0.5
      Transparency of the background image.
      0: Fully visible
      1: Fully transparent
    """
    data = self.getImageData()
    
    if np.all(data) == None:
      raise ValueError("No image loaded")

    if sigma is not None:
      data = gaussian_filter(data, sigma = sigma)

    norm = simple_norm(data, norm_type)

    fig = plt.figure()
    ax = WCSAxes(fig, [0, 0, 1, 1], wcs=self.getWCS())
    fig.add_axes(ax)

    if overlay:
      plt.imshow(self.getImageData(), cmap = base_cmap, alpha = alpha, 
        vmin=np.nanpercentile(self.getImageData(), 5), 
        vmax=np.nanpercentile(self.getImageData(), 99)
      )

    contour = plt.contour(data, levels = levels, cmap = cmap,
      vmin=np.nanpercentile(self.getImageData(), 5), 
      vmax=np.nanpercentile(self.getImageData(), 99)
    )
    cbar = plt.colorbar(contour)
    cbar.set_label("Intensity")

    plt.title("Contour plot")
    plt.xlabel("X-Axis")
    plt.ylabel("Y-Axis") 
    plt.show()
  
  def plotHist(self, nbins=200):
    # flatten means: we put our 2d array in a 1d array
    plt.hist(self.getImageData().flatten(), nbins)

    plt.xlabel('Pixel Content')
    plt.ylabel('Number of Pixels')
    plt.yscale('log')