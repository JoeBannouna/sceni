from astropy.io import fits 
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import proj_plane_pixel_scales, skycoord_to_pixel
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astroquery.vizier import Vizier

#Takes in parameters set by user and outputs stars within image corresponding to given star catalogue
class StarsFilter:

    def __init__ (self):
        self.catalogue_id = "I/239"
        self.mag_min = None
        self.mag_max = None
        self.period_min = None
        self.period_max = None


    #Sets star catalogue
    def set_catalogue(self, catalogue_id):
        self.catalogue_id = catalogue_id        
    
    #Gets star catalogue
    def get_catalogue(self):
        return self.catalogue_id
    
    #Defines region in right ascension and declination
    def set_region(self, ra_min, ra_max, dec_min, dec_max):
        self.ra_min = ra_min
        self.ra_max = ra_max
        self.dec_min = dec_min
        self.dec_max = dec_max
    
    #Gets region previously defined
    def get_region(self):
        return self.ra_min, self.ra_max, self.dec_min, self.dec_max
    
    #Sets limits for apparent magnitudes of stars found
    def set_mag_limit(self, mag_min, mag_max):
        self.mag_min = mag_min
        self.mag_max = mag_max
        return
    
    #Gets magnitude limits
    def get_mag_limit(self):
        return self.app_mag_min, self.app_mag_max
    
    #Sets limits for periods of stars found
    def set_periodicity_limit(self, period_min, period_max):
        self.period_min = period_min
        self.period_max = period_max
        
    #Gets period limits
    def get_periodicity(self):
        return self.period_min, self.period_max

    #Generates a list of stars which are within the given image with respect to the parameters given by the user previously
    #CURRENTLY ONLY WORKS WITH HIPPARCUS CATALOGUE
    def set_visible_stars(self, image):
        # Save header info
        wcs = WCS(image.header)
        
        # Obtain catalog using astroquery utilising the catalog ID
        Vizier.ROW_LIMIT = -1
        catalog_data = Vizier.get_catalogs(self.catalogue_id)[0]
        catalog_df = catalog_data.to_pandas()

        #Create a copy of the catalog DF where the only stars remaining are those that are within the bounds set by the image
        catalog_df = catalog_df[(catalog_df['_RA.icrs'] >= self.ra_min) &
                                (catalog_df['_RA.icrs'] <= self.ra_max) &
                                (catalog_df['_DE.icrs'] >= self.dec_min) &
                                (catalog_df['_DE.icrs'] <= self.dec_max)]

        #Save the coordinates of the filtered catalog and their unit (degrees) within a SkyCoord object
        star_coords = SkyCoord(ra = catalog_df['_RA.icrs'], dec = catalog_df['_DE.icrs'], unit = 'deg')
        
        #create 2 pixel arrays with pixels corresponding to star positions on the image
        x_pixels, y_pixels = skycoord_to_pixel(star_coords, wcs)

        #Add columns of xPixels and yPixels onto catalog
        catalog_df['x_pixels'] = x_pixels
        catalog_df['y_pixels'] = y_pixels

        catalog_df = catalog_df[(catalog_df['x_pixels'] >= 0) &
                                (catalog_df['x_pixels'] <= image.shape[1]) &
                                (catalog_df['y_pixels'] >= 0) &
                                (catalog_df['y_pixels'] <= image.shape[0])]

        #Check if apparent magnitude was set by user, and if so further filter stars according to the limits given
        if (self.mag_min != None):
                catalog_df = catalog_df[(catalog_df['Vmag'] >= self.mag_min)]
        if (self.mag_max != None):
                catalog_df = catalog_df[(catalog_df['Vmag'] <= self.mag_max)]

        # #Determine mathematically if stars are visible on the image by comparing surrounding pixel
        # for i in range(0, catalog_df.shape[0]):
        #     star_pixel_counts = image[catalog_df.iloc[catalog_df.columns.get_loc('x_pixels'), i]][catalog_df.iloc[catalog_df.columns.get_loc('y_pixels'), i]]
        #     print(star_pixel_counts)
        
        self.visible_stars = catalog_df
    
    #Returns the list of stars
    def get_visible_stars(self):
        return self.visible_stars