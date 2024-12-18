import os
import warnings
import numpy as np
import pandas as pd
from astropy.io import fits 
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier

from astrosceni.image import Image

#Simplifies warning message
warnings.simplefilter('always', UserWarning)
warnings.formatwarning = lambda message, *args: f"{message}\n"

#Takes in parameters set by user and outputs stars within image corresponding to given star catalogue
class StarsFilter:

    #Initializes class
    def __init__ (self):
        self.original_catalog_df = None
        self.visible_stars_df = None
        self.mag_min = None
        self.mag_max = None
        self.custom_region = False

    #Sets star catalogue
    #Uses hipparcos catalogue by default, can set catalogue id, and if custom catalogie user must input RA, Dec and apparent magnitude column names
    def setCatalogue(self, catalogue_id = "I/239/hip_main", ra_col_name = "_RA.icrs", dec_col_name = "_DE.icrs", app_mag_col_name = "Vmag"):

        #Check if given values are strings
        if (not isinstance(catalogue_id, str)) or (not isinstance(ra_col_name, str)) or (not isinstance(dec_col_name, str)) or (not isinstance(app_mag_col_name, str)):
            raise TypeError("Value must be a string.")

        #Save column names and save a copy of the catalog obtained
        self.ra_column_name = ra_col_name
        self.dec_column_name = dec_col_name
        self.app_mag_column_name = app_mag_col_name

        #Designate a file path
        file_name = catalogue_id.replace("/", "_") + ".csv"
        target_folder = "data"
        file_path = f"{target_folder}/{file_name}"

        if (os.path.isfile(file_path)):
            #Previous file generated was found, extracting dataframe
            print("Previous saved catalog file found.")
            catalog_df = pd.read_csv(file_path)
        else:
            #Previous file generated wasn't found, saving a new file with dataframe
            print("Previous saved catalog file not found, obtaining and saving new copy.")

            #Obtain catalog using astroquery utilising the catalog ID
            Vizier.ROW_LIMIT = -1
            catalog_data = Vizier.get_catalogs(catalogue_id)[0]
            catalog_df = catalog_data.to_pandas()
            catalog_df.to_csv(file_path, index = False)

        #Save catalog dataframe to class
        self.original_catalog_df = catalog_df.copy()
    
    #Returns original star catalog
    def getCatalogue(self):
        return self.original_catalog_df
    
    #Defines region in right ascension and declination
    def setRegion(self, ra_min, ra_max, dec_min, dec_max):

        #If region arguments are not floats, returns type error
        if (not isinstance(ra_min, float)) or (not isinstance(ra_max, float)) or (not isinstance(dec_min, float)) or (not isinstance(dec_max, float)):
            raise TypeError("Value must be a float.")

        #Saves region parameters, class remembers it now has a custom region
        self.custom_region = True
        self.ra_min = ra_min
        self.ra_max = ra_max
        self.dec_min = dec_min
        self.dec_max = dec_max
    
    #Returns region previously defined
    def getRegion(self):
        return self.ra_min, self.ra_max, self.dec_min, self.dec_max
    
    #Sets limits for apparent magnitudes of stars within image
    def setMagLimit(self, mag_min, mag_max):
        self.mag_min = mag_min
        self.mag_max = mag_max
    
    #Returns apparent magnitude limits
    def getMagLimit(self):
        return self.app_mag_min, self.app_mag_max
    
    #Sets limits for periods of stars found
    def setPeriodicityLimit(self, period_min, period_max):
        self.period_min = period_min
        self.period_max = period_max
        
    #Return period limits
    def getPeriodicityLimit(self):
        return self.period_min, self.period_max

    #Generates a list of stars which are within the given image with respect to the parameters given by the user previously
    def setVisibleStars(self, image):
        #Checks if image passed is an image object
        if (not isinstance(image, Image)):
            raise TypeError("Parameter must be an image object")
        
        #Checks that hdu from image isn't none
        hdu = image.getImageData(original = False)
        if (hdu is None):
            raise ValueError("hdu has datatype None")

        #If no catalogue currently loaded into program, use default catalogue, hipparcos
        if (self.original_catalog_df is None):
            warnings.warn("WARNING: No catalogue set, will use default star catalogue 'Hipparcos'")
            self.setCatalogue()

        #Create copy of catalog to alter
        catalog_df = self.original_catalog_df.copy()

        #If user never specified a custom region, will analyse entire image passed to it
        if self.custom_region == False:
            #Define the corners of the image
            corners = [
                (0, 0),                         # Bottom-left
                (hdu.shape[1], 0),              # Bottom-right
                (0, hdu.shape[0]),              # Top-left
                (hdu.shape[1], hdu.shape[0])    # Top-right
            ]

            #Decide the Ra and Dec range of the image
            raCorners, decCorners = (image.getWCS()).pixel_to_world_values(*zip(*corners))
            self.ra_min = min(raCorners)
            self.ra_max = max(raCorners)
            self.dec_min = min(decCorners)
            self.dec_max = max(decCorners)

        #Create a copy of the catalog DF where the only stars remaining are those that are within the bounds set by the image
        catalog_df = catalog_df[(catalog_df[self.ra_column_name] >= self.ra_min) &
                                (catalog_df[self.ra_column_name] <= self.ra_max) &
                                (catalog_df[self.dec_column_name] >= self.dec_min) &
                                (catalog_df[self.dec_column_name] <= self.dec_max)]

        #Save the coordinates of the filtered catalog and their unit (degrees) within a SkyCoord object
        star_coords = SkyCoord(ra = catalog_df[self.ra_column_name], dec = catalog_df[self.dec_column_name], unit = 'deg')
        
        #create 2 pixel arrays with pixels corresponding to star positions on the image
        x_pixels, y_pixels = skycoord_to_pixel(star_coords, image.getWCS())

        #Add columns of xPixels and yPixels onto catalog
        catalog_df['x_pixels'] = x_pixels
        catalog_df['y_pixels'] = y_pixels

        catalog_df = catalog_df[(catalog_df['x_pixels'] >= 0) &
                                (catalog_df['x_pixels'] <= hdu.shape[1]) &
                                (catalog_df['y_pixels'] >= 0) &
                                (catalog_df['y_pixels'] <= hdu.shape[0])]
        
        #Check if apparent magnitude was set by user, and if so further filter stars according to the limits given
        if (self.mag_min != None):
                catalog_df = catalog_df[(catalog_df[self.app_mag_column_name] >= self.mag_min)]
        if (self.mag_max != None):
                catalog_df = catalog_df[(catalog_df[self.app_mag_column_name] <= self.mag_max)]


        self.visible_stars_df = catalog_df
    
    #Returns the list of stars
    def getVisibleStars(self):
        return self.visible_stars_df