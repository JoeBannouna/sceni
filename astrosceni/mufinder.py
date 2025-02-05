import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astrosceni.image import Image

class MuFinder():

    def __init__(self, narrow_band_image, broad_band_image, mu_range=None, mu_resolution=0.1):
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
        self.mu_resolution = mu_resolution
        self.is_paramaters_changed = True   # Flag that the parameters have changes

    def setMuRange(self, start, end):
        """
        Sets the range of scale factors manually.
        """
        interval = self.mu_resolution

        assert start <= end, f"The minimum value, {start} is not less than the maximum value, {end}!"
        assert (end - start) > interval, f"The interval {interval} is not less than the range, {end - start}!"

        divisions = int((end - start)/interval)
        self.mu_linspace = np.linspace(start, end, divisions)
        self.is_paramaters_changed = True   # Flag that the parameters have changes

    def setMuRangeAuto(self):
        """
        Sets the range of scale factors automatically.
        """

        # INSERT CODE TO FIND THE BEST RANGE HERE

        self.mu_linspace = np.linspace(start, end, int((end - start)/self.mu_resolution))
        self.is_paramaters_changed = True   # Flag that the parameters have changes

    def getSkewnessVals():
        if (self.is_paramaters_changed):
            self.skewness_vals = []
            for mu in scale_factor_linspace:
                image_data = self.narrow_band_image.getImageData() - mu*self.broad_band_image.getImageData()
                data = np.array(image_data.flatten())
                skews.append(np.sum(((data - np.mean(data)) / np.std(data, ddof=1)) ** 3) / (len(data) - 1))

        return self.skewness_vals

    def plotSkewnessVals():
        plt.plot(self.factors, self.getSkewnessVals())

        # Set axis labels, anything in between $^\$ will be displayed as a mathematical symbol
        plt.xlabel('$^\mu$')
        plt.ylabel('s($^\mu$)')
        

    # def (self):
        


    # def plotPixelDistribution(self, scale_factor = 0.1, nbins = 301, lower_lim = -45000.0, upper_lim = 45000.0):
    #     """
    #     Plots a histogram of the pixel distribution for a specific scale factor.

    #     Parameters
    #     ----------
    #     scale_factor: float
    #         The scale factor to consider.
    #         Optional input.
    #         Default value is 0.1.
        
    #     nbins: integer
    #         Number of bins for the histogram.
    #         Optional input.
    #         Default value is 301.
        
    #     lower_lim: int, float
    #         Lower limit for the range of the x-axis for the histogram.
        
    #     upper_lim: int, float
    #         Upper limit for the range of the x-axis for the histogram.
    #     """

    #     self.residual_image = Subtractor.getSubtractedImage(self.NB_arr, self.BB_arr, scale_factor)

    #     # Flatten the array for histogram plotting
    #     flat_residual = self.residual_image.getImageData().flatten()

    #     # Plot the histogram
    #     plt.hist(flat_residual, bins = nbins, range = (lower_lim, upper_lim))

    #     plt.text(0.5 * upper_lim, 0.9 * plt.ylim()[1], f"mu = {scale_factor}", ha='left', va='top')
    #     plt.text(0.5 * upper_lim, 0.4 * plt.ylim()[1], f"skew = {Subtractor.skewnessCalculator(self.residual_image.getImageData()): .2f}", ha='left', va='top')
    #     plt.yscale('log')

    #     # Label axes and display the plot
    #     plt.xlabel("Pixel Value")
    #     plt.ylabel("Frequency")
    #     plt.title("Pixel Distribution")
    #     plt.show()


    # def plotSkewVsMu(self):
    #     plt.plot(self.factors, self.skewness_vals)

    #     # Set axis labels, anything in between $^\$ will be displayed as a mathematical symbol
    #     plt.xlabel('$^\mu$')
    #     plt.ylabel('s($^\mu$)')

    
    # def getResultImage(self):
    #     # Turn Residual array back into a image class.
    #     pass