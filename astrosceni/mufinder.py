import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astrosceni.image import Image
from scipy.interpolate import make_splrep, PPoly
from scipy.stats import skew

class MuFinder:
    def __init__(self, narrow_band_image, broad_band_image, mu_range=None, mu_resolution=0.1):
        """
        Instantiates an instance of the MuFinder class.

        Parameters
        ----------
        narrow_band_image: image file
            The narrow band image file.

        broad_band_image: image file
            The broad band image file.

        mu_range: tuple, list, array
            (min, max): The minimum and maximum values to consider when finding the optimum mu.
            Optional
            Default: Applies the setMuRangeAuto function.
        
        mu_resolution: float
            Optional 
            Default: 0.1
            Increment between each mu value to consider within the specified range.
        """
        # Class instance state
        self.narrow_band_image = narrow_band_image
        self.broad_band_image = broad_band_image
        self.mu_range = mu_range
        self.mu_resolution = mu_resolution

        self.is_paramaters_changed = True
        self.mu_linspace = None
        self.skewness_vals = None
        self.optimal_mu = None
        self.result_image = None

        if (mu_range == None): 
            self.setMuRangeAuto()
        else:
            divisions = int((self.mu_range[1] - self.mu_range[0])/self.mu_resolution)
            self.mu_linspace = np.linspace(self.mu_range[0], self.mu_range[1], divisions)


    def setMuResolution(self, mu_resolution):
        """
        Sets the increment between each mu value to consider within the specified range.
        """
        self.mu_resolution = mu_resolution
        self.is_paramaters_changed = True   # Flag that the parameters have changes


    def setMuRange(self, start, end):
        """
        Sets the range of scale factors manually.

        Parameters
        ----------
        start: float
            Minimum mu value.

        end: float
            Maximum mu value.
        """
        interval = self.mu_resolution

        # Restrict end >= start
        assert start <= end, f"The minimum value, {start} is not less than the maximum value, {end}!"

        # Restrict end - start to be > interval
        assert (end - start) > interval, f"The interval {interval} is not less than the range, {end - start}!"

        divisions = int((end - start)/interval)
        self.mu_linspace = np.linspace(start, end, divisions)
        self.is_paramaters_changed = True   # Flag that the parameters have changes


    def setMuRangeAuto(self):
        """
        Sets the range of scale factors automatically.

        1. Gets a copy of the narrow and broad band image data.
        2. Divides each pixel value of the narrowband data by the corresponding broadband data.
        3. Sets the median value of these ratios to be the median mu.
        4. Chooses the larger between 0.1 and 0.5*(median mu) and sets this to be the min mu.
        5. Sets the max mu to be 2*(median mu)
        """

        narrow_data = self.narrow_band_image.getImageData()
        broad_data = self.broad_band_image.getImageData()

        # Safeguard to make sure not dividing by 0
        epsilon = 1e-6
        valid_pixels = broad_data > epsilon

        # Determining the rough ratio of the two images, obtaining median mu
        pixel_ratios = narrow_data[valid_pixels] / (broad_data[valid_pixels] + epsilon)
        median_mu = np.nanmedian(pixel_ratios)

        # Based on rough median_mu, choosing large enough range
        start = max(0.1, 0.5 * median_mu)
        end = (2 * median_mu)

        self.mu_range = (start, end)
        self.mu_linspace = np.linspace(start, end, int((end - start)/self.mu_resolution))
        self.is_paramaters_changed = True   # Flag that the parameters have changes


    def getSkewnessVals(self):
        """
        Get the skewness values for each mu value considered.

        Returns
        -------
        skewness_vals: array
            An array of the skewness values for each corresponding mu value.
        """
        if (self.is_paramaters_changed):
            self.skewness_vals = []
            narrow_data = self.narrow_band_image.getImageData().flatten()
            broad_data = self.broad_band_image.getImageData().flatten()
            
            bad_data = np.isnan(narrow_data) | np.isnan(broad_data)
            narrow_data = narrow_data[~bad_data]
            # narrow_data = narrow_data[~np.isnan(broad_data)]
            broad_data = broad_data[~bad_data]
            # broad_data = broad_data[~np.isnan(narrow_data)]

            for mu in self.mu_linspace:
                data = narrow_data - mu * broad_data  # Element-wise operation
                skew_val = skew(data, bias=False)
                self.skewness_vals.append(skew_val)
        self.is_paramaters_changed = False
        return self.skewness_vals


    def plotSkewnessVals(self, h_line = True, v_line = True):
        """
        Plots the skewness vs the corresponding mu scaling factors

        Parameters:
        -------
        h_line: boolean
            optional
            default: True
            Draws a horizontal line where skewness = 0
        
        v_line: boolean
            optional
            default: True
            Draws a vertical line at the Mu closest to where skewness = 0

        """
        skewnessVals = self.getSkewnessVals()
        plt.xlabel('$^\mu$')
        plt.ylabel('s($^\mu$)')
        plt.plot(self.mu_linspace, skewnessVals)

        if h_line == True:
            plt.axhline(y = 0)

        if v_line == True:
            plt.axvline(x = self.getCorrespondingMu())

    def getCorrespondingMu(self, skewness=0):
        """
        Returns the corresponding mu for a given skew

        Parameters:
        -------
        skewness: float
            optional
            default: 0
            The mu that will be returned is the mu in the mu_linspace that best matches to this skewness
        """
        skewness_vals = self.getSkewnessVals()
        closest_idx = np.argmin(np.abs(np.array(skewness_vals)-skewness))

        return self.mu_linspace[closest_idx]


    def getOptimalMus(self):
        """
        Determines the optimal scaling factor (mu).

        Returns
        -------
        roots: float
            The optimal scaling factor (mu).
        """

        # Obtains possible optimal mus.
        # The smaller of the possible values is the true optimal mu.

        skews = self.getSkewnessVals() # This is a numpy array
        roots = PPoly.from_spline(make_splrep(self.mu_linspace, skews).derivative()).roots() # Roots are the optimal mu. Use spline to find derivative to find optimal mu.
        # PPoly.from_spline constructs a piecewise polynomial from a spline
        # numpy.roots(p) outputs a numpy array with the roots of a polynomial with coefficients p.
        # e.g. roots(1, 3, 2) outputs array([-1,, -2])
        # Access the coefficients of PPoly.from_spline(make_splrep(self.mu_linspace, skews).derivative()) using .c e.g.:
        #   PPoly.from_spline(make_splrep(self.mu_linspace, skews).derivative()).c
        
        roots = roots[roots < self.mu_range[1]]
        roots = roots[roots > self.mu_range[0]]

        #print(roots2ndDerivative)
        return roots


    def getResultImages(self):
        """
        Gets the result of performing NB - mu*BB where NB is the narrowband image, BB is the broadband image and mu is the optimal scale factor.

        Returns
        -------
        images: array
            Resultant images corresponding to all possible optimal mus determined from the getOptimalMus function.
            Since the mu with the lowest value among the possible optimal mus is the true optimal scale factor, the first image is that which is the result of subtracting using the true optimal mu.
            Access different images using "image[i]" where i is the index number.
        """
        # Gets the resultant image for each corresponding mu
        # Lowest mu corresponds to first image in array.
        # Image with index 0 is the resultant image using the true optimal mu.
        # can get image using images[i] where i is the index.
        optimal_mus = self.getOptimalMus()
        images = []
        for mu in optimal_mus:
            images.append(Image.subtract(self.narrow_band_image, self.broad_band_image, mu))
        return images