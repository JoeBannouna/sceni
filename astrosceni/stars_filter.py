import copy
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
from scipy.optimize import curve_fit
from astrosceni.image import Image

#Simplifies warning message
warnings.simplefilter('always', UserWarning)
warnings.formatwarning = lambda message, *args: f"{message}\n"

#Takes in parameters set by user and outputs stars within image corresponding to given star catalogue
class StarsFilter:

    #Initializes class
    def __init__ (self, data_directory_path = 'data'):
        """
        Instantiates an instance of the StarsFilter class.
        
        Parameters
        ----------
        data_directory_path: string
            Optional
            Default: 'data' (This assumes the home directory is the 'SCENI' directory.)
            The path to the directory where the star catalogue files are situated in. 
            Note that this parameter only takes in the path to a directory not a specific file as all the files in said directory will be read.
        """
        self.original_catalog_df = None
        self.stars_in_region_df = None
        self.stars_visible_df = None
        self.mag_min = None
        self.mag_max = None
        self.custom_region = False
        self.catalogue_set = False

        self.data_directory_path = data_directory_path 


    def setCatalogue(self, download_catalogue = True, catalogue_id = "I/239/hip_main", ra_col_name = "_RA.icrs", dec_col_name = "_DE.icrs", app_mag_col_name = "Vmag"):
        """
        Sets star catalogue desired, and by default saves copy of important star properties (Ra, Dec and Vmag of stars).
        Uses hipparcos catalogue by default.
        Do note that different star catalogues use different column names thus, do not always assume that the default is correct.

        Parameters
        ----------
        download_catalogue: boolean
            Optional
            Default: True
            Controls whether or not to download the star catalogue files.

        catalogue_id: string
            Optional
            Default: "I/239/hip_main"
            The star catalogue file to use. Only the file name is required, thus exclude the file type, e.g. ".csv".

        ra_col_name: string
            Optional
            Default: "_RA.icrs" (This is for the Hipparcos star catalogue)
            Column name for column with right ascension coordinates.

        dec_col_name: string
            Optional
            Default: "_DE.icrs" (This is for the Hipparcos star catalogue)
            Column name for the column with the declination coordinates.

        app_mag_col_name: string
            Optional
            Default: "Vmag" (This is for the Hipparcos star catalogue)
            Column name for the column with the apparent magnitudes.
        """
        #connects to data base online, finds the catalogue based on ID.
        # Saves col names you want only
        # put in col names and it will only return the columns
        # Different catalogues have different column names 

        #Check if given values are strings
        # Use the column names as shown on the online databse.
        if (not isinstance(catalogue_id, str)) or (not isinstance(ra_col_name, str)) or (not isinstance(dec_col_name, str)) or (not isinstance(app_mag_col_name, str)):
            raise TypeError("Value must be a string.")

        #Save column names and save a copy of the catalog obtained
        self.ra_col_name = ra_col_name.replace("(", "_").replace(")", "_").replace("-", "_")
        self.dec_col_name = dec_col_name.replace("(", "_").replace(")", "_").replace("-", "_")
        self.app_mag_col_name = app_mag_col_name.replace("(", "_").replace(")", "_").replace("-", "_")

        #Designate a file path
        file_name = catalogue_id.replace("/", "_") + ".csv"
        target_folder = self.data_directory_path
        file_path = f"{target_folder}/{file_name}"

        if (os.path.isfile(file_path)):
            print("Previous saved catalog file found.")
            catalog_df = pd.read_csv(file_path)
        else:
            print("Previous saved catalog file not found")
            print("Downloading new copy")

            #Obtain catalog using astroquery utilising the catalog ID
            vizier = Vizier(columns = [ra_col_name, dec_col_name, app_mag_col_name])
            vizier.ROW_LIMIT = -1
            catalog_data = vizier.get_catalogs(catalogue_id)[0]
            catalog_df = catalog_data.to_pandas()
            print(catalog_df.columns)

            catalog_df = catalog_df.rename(columns={catalog_df.columns[0]: "RA", catalog_df.columns[1]: "DEC", catalog_df.columns[2]: "Vmag"})

            #Unless user doesnt specify, download catalogue to file
            if download_catalogue == True:
                catalog_df.to_csv(file_path, index = False)
                print("Saving copy of catalogue from online into .csv file")
            else:
                print("Obtaining copy of catalogue from online")

        self.catalogue_set = True # Prevents resetting the catalogue at the end of the function.

        #Save catalog dataframe to class
        self.original_catalog_df = catalog_df.copy()
    

    #Returns original star catalog
    def getCatalogue(self):
        """
        Returns the original star catalog.
        
        Returns
        -------
        original_catalog_df: Pandas Dataframe
            The original star catalog."""
        return self.original_catalog_df
    

    #Defines region in right ascension and declination
    def setRegion(self, ra_min, ra_max, dec_min, dec_max):
        """
        Defines region in right ascension and declination.

        Parameters
        ----------
        ra_min: float
            Minimum right ascension coordinate for the region.

        ra_max: float
            Maximum right ascension coordinate for the region.

        dec_min: float
            Minimum declination coordinate for the region.

        dec_max: float
            Maximum declination coordinate for the region.
        """
        # If region arguments are not floats, returns type error
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
        """
        Returns the region defined in the setRegion function.

        Returns
        -------
        ra_min: float
            Minimum right ascension coordinate for the region.

        ra_max: float
            Maximum right ascension coordinate for the region.

        dec_min: float
            Minimum declination coordinate for the region.

        dec_max: float
            Maximum declination coordinate for the region.
        """
        return self.ra_min, self.ra_max, self.dec_min, self.dec_max
    
    #Sets limits for apparent magnitudes of stars within image
    def setMagLimit(self, mag_min, mag_max):
        """
        Sets the limits for the apparent magnitudes of stars within the image.

        Parameters
        ----------
        mag_min: float
            The minimum apparent magnitude.

        mag_max: float
            The maximum apparent magnitude.
        """
        self.mag_min = mag_min
        self.mag_max = mag_max
    

    #Returns apparent magnitude limits
    def getMagLimit(self):
        """
        Returns the apparent magnitude limits as defined in the setMagLimit function.

        Returns
        ----------
        mag_min: float
            The minimum apparent magnitude.

        mag_max: float
            The maximum apparent magnitude.
        """
        return self.mag_min, self.mag_max
    

    #Sets limits for periods of stars found
    def setPeriodicityLimit(self, period_min, period_max):
        """
        Sets the limits for the periods of stars found.

        Parameters
        ----------
        period_min: float
            Minimum period.

        period_max: float
            Maximum period.
        """
        self.period_min = period_min
        self.period_max = period_max
        

    #Return period limits
    def getPeriodicityLimit(self):
        """
        Returns the limits for the periods of stars as defined in the setPeriodicityLimits function.

        Returns
        -------
        period_min: float
            Minimum period.

        period_max: float
            Maximum period.
        """
        return self.period_min, self.period_max


    def setStarsInRegion(self, image):
        """
        Generates dataframe of stars which are within given image using header information
        Dataframe Layout:
        Star Identifier #       RA      DEC     Vmag        x_pixels        y_pixels
        ...                     ...     ...     ...         ...             ...     

        Parameter
        ---------
        image: numpy array
            The numpy array of the data for the image. This should be an instance of the Image class - an image object.
        """

        # Checks if image passed is an image object
        if (not isinstance(image, Image)):
            raise TypeError("Parameter must be an image object")
        
        # Checks that data from image isn't none
        data = image.getImageData(original = False)
        if (data is None):
            raise ValueError("data has datatype None")

        # If no catalogue currently loaded into program, use default catalogue, hipparcos
        if (self.original_catalog_df is None):
            warnings.warn("WARNING: No catalogue set, will use default star catalogue 'Hipparcos'")
            self.setCatalogue()

        # Create copy of catalog to alter
        catalog_df = self.original_catalog_df.copy()

        # If user never specified a custom region, will analyse entire image passed to it
        if self.custom_region == False:
            #Define the corners of the image
            corners = [
                (0, 0),                         # Bottom-left
                (data.shape[1], 0),              # Bottom-right
                (0, data.shape[0]),              # Top-left
                (data.shape[1], data.shape[0])    # Top-right
            ]

            # Decide the Ra and Dec range of the image
            raCorners, decCorners = (image.getWCS()).pixel_to_world_values(*zip(*corners))
            self.ra_min = min(raCorners)
            self.ra_max = max(raCorners)
            self.dec_min = min(decCorners)
            self.dec_max = max(decCorners)

        # Create a copy of the catalog DF where the only stars remaining are those that are within the bounds set by the image
        catalog_df = catalog_df[(catalog_df["RA"] >= self.ra_min) &
                                (catalog_df["RA"] <= self.ra_max) &
                                (catalog_df["DEC"] >= self.dec_min) &
                                (catalog_df["DEC"] <= self.dec_max)]

        # Save the coordinates of the filtered catalog and their unit (degrees) within a SkyCoord object
        star_coords = SkyCoord(ra = catalog_df["RA"], dec = catalog_df["DEC"], unit = 'deg')

        # create 2 pixel arrays with pixels corresponding to star positions on the image
        x_pixels, y_pixels = skycoord_to_pixel(star_coords, image.getWCS())

        # Add columns of xPixels and yPixels onto catalog
        catalog_df['x_pixels'] = x_pixels
        catalog_df['y_pixels'] = y_pixels

        catalog_df = catalog_df[(catalog_df['x_pixels'] >= 0)]
        catalog_df = catalog_df[(catalog_df['x_pixels'] <= data.shape[1])]
        catalog_df = catalog_df[(catalog_df['y_pixels'] >= 0)]
        catalog_df = catalog_df[(catalog_df['y_pixels'] <= data.shape[0])]

        #Check if apparent magnitude was set by user, and if so further filter stars according to the limits given
        if (self.mag_min != None):
                catalog_df = catalog_df[(catalog_df["Vmag"] >= self.mag_min)]
        if (self.mag_max != None):
                catalog_df = catalog_df[(catalog_df["Vmag"] <= self.mag_max)]

        self.stars_in_region_df = catalog_df


    def extractStarRegion(self, star_index, image, x_y_width):
        """
        Extracts an n x n grid around the centre of a star for analysis.

        Parameters
        ----------
        star_index: integer
            Index of a star in the star catalogue file.

        image: numpy array
            The numpy array of the data of the image. This is an instance of the Image class - an image object.

        x_y_width: integer 
            The halfwidth of the square region to be extracted.

        Returns
        -------
        padded_region: numpy array
            The region which extends beyond the image boundaries. This region will be padded with NaN values.
        """
        star = self.stars_in_region_df.iloc[star_index]
        x, y = int(star['x_pixels']), int(star['y_pixels'])
        data = image.getImageData()

        # Define the bounds of the region
        x_min = max(x - x_y_width, 0)
        x_max = min(x + x_y_width, data.shape[1])
        y_min = max(y - x_y_width, 0)
        y_max = min(y + x_y_width, data.shape[0])

        # Extract the region
        region = data[y_min:y_max, x_min:x_max]

        # Pad the region with NaN values if it extends beyond the image boundaries
        padded_region = np.full((2 * x_y_width, 2 * x_y_width), np.nan)
        padded_region[(y_min - y + x_y_width):(y_max - y + x_y_width), (x_min - x + x_y_width):(x_max - x + x_y_width)] = region

        return padded_region


    def plotHistOfStar(self, star_index, image1, image2, print_brightest_pixels = False):
        """
        Plots a histogram of the star region, with both image 1 (NB) and image 2 (BB).

        Parameters
        ----------
        star_index: integer
            Index of a star in the star catalogue file.

        image1: numpy array
            The data for the narrowband image. This should be an instance of the image class - an image object.

        image2: numpy array
            The data for the broadband image. This should be an instance of the image class - an image object.

        print_brightest_pixels: boolean
            Controls whether or not to print a line each indicating the brightest pixel within the specified range for the narrowband and broadband images.
        """

        region1 = self.extractStarRegion(star_index, image1, 5)
        region2 = self.extractStarRegion(star_index, image2, 5)

        if print_brightest_pixels == True:
            print("IMAGE 1, Brightest pixel within range: ", np.nanmax(region1))
            print("IMAGE 2, Brightest pixel within range: ", np.nanmax(region2))
    
        plt.figure()
        plt.hist(region1.ravel(), bins=30, alpha = 0.7, color = 'red', label = 'NB Pixel Values')
        plt.hist(region2.ravel(), bins=30, alpha = 0.7, color = 'blue', label = 'BB Pixel Values')
        plt.legend()
        plt.xlabel("Pixel Value")
        plt.ylabel("Frequency")
        plt.title(f"Histogram of pixel values around star with index {star_index}")


    def determineVisibilityOfIndividual(self, star_index, image, region_size=10, wiggleRoom = 2, threshold_amp = 500, reduced_chi_squared_bound = 2, print_results = False):
        """
        Determines whether a star is visible based on whether a gaussian can be fitted to it, and if the gaussian has a high enough amplitude and low enough reduced chi squared.

        Parameters
        ----------
        star_index: integer
            Index of a star in the star catalogue file.

        image: numpy array
            The image data. This should be an instance of the image class - an image object.

        region_size: integer
            Optional
            Default: 10
            Size of the region to consider.

        wiggleRoom: integer
            Optional
            Default: 2
            The size of the region to check for the center of a star.
            This is in case the position of a star in the star catalogue is slightly off.I
            If catalogue is slightly inaccurate, the wiggle room checks if there's any oyher brightr pixels around it and shifts the point considered to be the center of the star.

        threshold_amp: integer
            Optional
            Default: 500
            Sets the amplitude condition to determine whether a star is visible.
            Amplitude is the apparent magnitude at the center of a star according to the Gaussian fit.
            Anything lower than 500 is definitely not visible.

        reduced_chi_squared_bound: integer
            Optional
            Default: 2
            Sets the chi squared condition to determine whether a star is visible.
            The lower the chi squared, the better the Gaussian fit.
            Through empirical testing, Gaussian fits usually have a chi square of 1 or less. Hence, a chi squared >2 is usually a really bad fit.

        print_results: boolean
            Optional
            Default: False
            Controls whether or not to print the result stating whether or not the star is visible.
        """

        region = self.extractStarRegion(star_index, image, region_size)

        # Check if the region is empty
        if np.all(np.isnan(region)):
            if print_results == True:
                print(f"Star index {star_index} has an empty region.")
            return False

        # Create grid of x and y coordinates which correspond to padded region pixels
        x = np.arange(region.shape[1])
        y = np.arange(region.shape[0])
        x, y = np.meshgrid(x, y)

        # Retrieve known star position (from cataloge) in full image
        star = self.stars_in_region_df.iloc[star_index]
        x_full = int(star['x_pixels'])
        y_full = int(star['y_pixels'])

        #Compute the extraction bounds (same as in extractStarRegion)
        x_min = max(x_full - region_size, 0)
        y_min = max(y_full - region_size, 0)

        # Initial guess for gaussian center as local maximum
        # Assume brightest point at middle of a star (best place to fit Guassian)
        x0_initial = x_full - x_min
        y0_initial = y_full - y_min

        # Define full 2D Gaussian function
        def gaussian_2d(coords, x0, y0, amplitude, sigma, offset):
            x, y = coords
            return amplitude * np.exp(-(((x-x0)**2)/(2*(sigma**2)) + ((y-y0)**2)/(2*(sigma**2)))) + offset

        # Initial guess for free parameters (amplitude, sigma, offset)
        initial_guess = [x0_initial, y0_initial, np.nanmax(region) - np.nanmedian(region), 2, np.nanmedian(region)]

        # Create a mask to ignore NaN valuyes in the region
        mask = ~np.isnan(region)

        #Bounds
        # x0 and y0 are bounded to be within 2 pixels of the initial guess
        bounds = (
        [x0_initial - wiggleRoom, y0_initial - wiggleRoom, 0, 0, -np.inf],
        [x0_initial + wiggleRoom, y0_initial + wiggleRoom, np.inf, np.inf, np.inf]
        )                  
        # changes few parameters of the gaussian for best fit
        # Gaussian found would have an amplitude
        # Computes a reduced chi square - used to determine how good a Gaussian fit is

        # Testing Gaussian
        try:
            popt, pcov = curve_fit(
                gaussian_2d,
                (x[mask], y[mask]),
                region[mask],
                p0=initial_guess,
                bounds=bounds
            )
        except RuntimeError:
            if print_results == True:
                print("Gaussian fitting failed")
            return False

        # Extract amplitude from fit parameters
        amplitude = popt[2]

        #Compute model values for entire region using fitted parameters
        model_values = gaussian_2d((x, y), *popt).reshape(region.shape)
        residuals = region - model_values

        #Calculate chi-squared using only valid datapoints
        chi_squared = np.sum((residuals[mask]/np.nanstd(region[mask])) ** 2)
        reduced_chi_squared = chi_squared / (mask.sum() - len(popt))

        # Check if star is bright enough and ensure fit is reliable
        if print_results == True:
            print("Star index: ", star_index, ". Amplitude: ", amplitude, ". Reduced Chi Squared: ", reduced_chi_squared)

        # Amplitude is the brightness of the middle pixel according to the Gaussian.
        # anything lower than 500, definitely not visible
        # THe lower the reduced chi square, the more accurate
        # Through empirical testing, usually most Gaussian fittings are around 1 or less, usually really bad Gaussian fitting will have chi squared >2
        if (amplitude > threshold_amp) and (reduced_chi_squared < reduced_chi_squared_bound):
            if print_results == True:
                print("Star is considered visible")
            return True # Star is visible
        
        if print_results == True:
            print(f"Star is not considered visible, amplitude is not greater than {threshold_amp} and reduced chi squared is not less than {reduced_chi_squared_bound}}")
        return False    # Star is not visible 


    def setVisibleStars(self, image, print_results = False, region_size = 10, wiggleRoom = 2, threshold_amp = 500, reduced_chi_squared_bound = 2):
        """
        Iterates through all stars in image and determines if they are visible, saves new, shortened, dataframe.

        Parameters
        ----------
        image: numpy array
            The image data. This should be an instance of the image class - an image object.

        print_results: boolean
            Optional
            Default: False
            Controls whether or not to print the result stating whether a star is visible.

        region_size: integer
            Optional
            Default: 10
            Sets the size of the region to consider.

        wiggleRoom: integer
            Optional
            Default: 2
            The size of the region to check for the center of a star.
            This is in case the position of a star in the star catalogue is slightly off.I
            If catalogue is slightly inaccurate, the wiggle room checks if there's any oyher brightr pixels around it and shifts the point considered to be the center of the star.

        threshold_amp: integer
            Optional
            Default: 500
            Sets the amplitude condition to determine whether a star is visible.
            Amplitude is the apparent magnitude at the center of a star according to the Gaussian fit.
            Anything lower than 500 is definitely not visible.

        reduced_chi_squared_bound: integer
            Optional
            Default: 2
            Sets the chi squared condition to determine whether a star is visible.
            The lower the chi squared, the better the Gaussian fit.
            Through empirical testing, Gaussian fits usually have a chi square of 1 or less. Hence, a chi squared >2 is usually a really bad fit.
        """
        visible_stars = []
        for i in range(len(self.stars_in_region_df)):
            if self.determineVisibilityOfIndividual(i, image, region_size, wiggleRoom, threshold_amp, reduced_chi_squared_bound, print_results):
                visible_stars.append(self.stars_in_region_df.iloc[i])
        self.stars_visible_df = pd.DataFrame(visible_stars)


    def removeVisibleStars(self, image, region_size = 10, annulus_width = 2):
        """
        Replaces visible stars with the estimated background level, size of which can be determined by user.

        Parameters
        ----------
        image: numpy array
            The image data. This should be an instance of the image class - an image object.
        
        region_size: integer
            Optional
            Default: 10
            Sets the size of the region to consider.

        annulus_width: integer
            Optional
            Default: 2
            Sets the width of the annulus around a star region to estimate the background brightness.

        Returns
        -------
        result: numpy array
            The new image data. This is an instance of the image class - an image object.
        """
        #Creates copy of image data
        new_image_data = image.getImageData().copy()

        #Global background median
        global_background_level = np.median(new_image_data)
        #Iterates through all visible stars in image and replaces them with background level
        for idx, star in self.stars_visible_df.iterrows():
            # Get stars pixel coordinates (already in dataframe)
            x_center, y_center = int(star["x_pixels"]), int(star["y_pixels"])

            #Define star region (a square centered on the star)
            x_min = max(x_center - region_size, 0)
            x_max = min(x_center + region_size, new_image_data.shape[1])
            y_min = max(y_center - region_size, 0)
            y_max = min(y_center + region_size, new_image_data.shape[0])

            #Define annulus region around star region for estimating the background
            #Expanding box by annulus_width on all sides
            ann_x_min = max(x_min - annulus_width, 0)
            ann_x_max = min(x_max + annulus_width, new_image_data.shape[1])
            ann_y_min = max(y_min - annulus_width, 0)
            ann_y_max = min(y_max + annulus_width, new_image_data.shape[0])

            #Create mask to ignore inner star region
            annulus = new_image_data[ann_y_min:ann_y_max, ann_x_min:ann_x_max]
            mask = np.ones_like(annulus, dtype=bool)

            #Mark star region as False
            inner_x_start = x_min - ann_x_min
            inner_x_end = inner_x_start + (x_max - x_min)
            inner_y_start = y_min - ann_y_min
            inner_y_end = inner_y_start + (y_max - y_min)
            mask[inner_y_start:inner_y_end, inner_x_start:inner_x_end] = False

            #Compute background median from annulus
            background_pixels = annulus[mask]
            if background_pixels.size > 0:
                background_level = np.median(background_pixels)
            else:
                # Fallback in case annulus is empty
                background_level = global_background_level

            # Replace star region with abckground level
            new_image_data[y_min:y_max, x_min:x_max] = background_level


        result = copy.deepcopy(image)
        result.setImageData(new_image_data)
        return result


    # Returns the list of stars in region
    def getStarsInRegion(self):
        """
        Returns the list of stars in a region
    
        Returns
        -------
        stars_in_region_df: pandas dataframe
        List of stars in a region defined in the extractStarRegion function.
        """

        return self.stars_in_region_df
    

    # Returns list of visible stars in region
    def getVisibleStars(self):
        """
        Returns a list of visible stars in a region
        
        Returns
        -------
        stars_visible_df: pandas dataframe
            List of visible stars in a region defined in the setVisibleStars function.
        """
        return self.stars_visible_df
    

    # Simple usage, pass the image and will return array of data with visible stars removed
    def filterStars(self, image, print_findings = False):
        """
        All encompassing function, acts as simple usage of class. Returns data of image with visible stars already removed, uses default parameters.

        Parameters
        ----------
        image: numpy array
            The image data. This should be an instance of the image class - an image object.

        print_findings: boolean
            Optional
            Default: False
            Controls whether or not to print the stars in the stipulated region and all the visible stars in the stipulated region.

        Returns
        -------
        A numpy array with all the visible stars replaced with the estimated background brightness. This is an instance of the image class - an image object.
        """    
        #If no catalogue specifically chosen, fall back to hipparcos
        if self.catalogue_set == False:
            self.setCatalogue(download_catalogue = True)
            self.catalogue_set = True

        # Sets dataframe of stars in region
        self.setStarsInRegion(image)

        # Sets dataframe of visible stars in region
        self.setVisibleStars(image)

        # Returns data of image with visible stars removed
        if print_findings == True:
            print("Stars in region:", len(self.stars_in_region_df))
            print("Visible stars:", len(self.stars_visible_df))

        return self.removeVisibleStars(image)