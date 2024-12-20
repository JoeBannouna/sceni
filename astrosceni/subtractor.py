import numpy as np

import pandas as pd

import matplotlib.pyplot as plt

from astrosceni.image import Image

class Subtractor():

    def __init__(self, NB, BB):
        """
        Parameters:
        -----------
        NB: image object
            The narrowband image object.
            
        BB: image object
            The broadband image object.
        """
  
        # convert NB image object to pandas dataframe
        self.NB_arr = np.array(NB.getImageData())


        # convert BB image object to pandas dataframe
        self.BB_arr = np.array(BB.getImageData())

   
    def setScaleFactorRange(self, min = 0.0, max = 2.0, interval = 0.01):
        """
        Sets the range of scale factors to determine the optimal scale factor.

        Parameters:
        ----------
        min: int, float
            Minimum scale factor value.
            Optional input.
            Default value is 0.

        max: int, float
            Maximum scale factor value.
            Optional input.
            Default value is 2.

        interval: float
            Interval between each scale factor to evaluate.
            Optional input.
            Default value is 0.01.
        """

        assert min <= max, f"The minimum value, {min} is not less than the maximum value, {max}!"
        assert (max - min) > interval, f"The interval {interval} is not less than the range, {max - min}!"

        divisions = int((max - min)/interval)

        self.factors = np.linspace(min, max, divisions)


    coord_types = ["coordinate system", "pixel value"]

    def setTestRegion(self, x_1, y_1, x_2, y_2, coord_type = "pixel value"):
        """
        This function is optional and is to be used if only a certain region in the image is to be considered to determine the optimal scale factor.

        Parameters
        ----------
        x_1, x_2: float
            Starting coordinates and ending coordinates.
            
        y_1, y_2: float
            Starting coordinates and ending coordinates.
        
        coord_type: string
            The coordinate type, allowed inputs are: "coordinate system" and "pixel value".
            With "coordinate system", the x_1, y_1, x_2, and y_2 parameters are as defined above.
            With "pixel value", the x_1, y_1, x_2, and y_2 parameters will be treated as pixel indices.
            Optional input.
            Default is "pixel value".
        """
        

        if coord_type == "coordinate system":
            pass # pass for now since we haven't figured out how to get ra and dec from a FITS file


        if coord_type == "pixel value":
            self.NB_arr = self.NB_arr[x_1: x_2, y_1: y_2]
            self.BB_arr = self.BB_arr[x_1: x_2, y_1: y_2]

    
    # define a function to calculate NB - mu*BB
    @staticmethod
    def getSubtractedImage(NB_image, BB_image, mu):
        """
        Calculates Narrowband - (scaling parameter)(Broadband).

        Parameters
        ----------
        NB_arr: numpy array 
            Narrowband numpy array.
        BB_arr: numpy array
            Broadband numpy array.
        
        mu: int, float
            The scale factor.
        
        Returns
        -------
        Residual: numpy array
            The residual numpy array"""

        if NB_image.getImageData().shape != BB_image.getImageData().shape: raise ValueError('Images must be the same size')
        result = Image.subtract(NB_image, BB_image, mu=mu)
        return result


    @staticmethod
    def skewnessCalculator(Residual_arr):
        """
        Calculates the skewness according to equation 1 from the paper (Sungryong Hong et al, 2014, Publications of the Astronomical Society of the Pacific).

        Parameters
        ----------
        Residual_arr: numpy array
            The numpy array of the result from executing the continuum_component_remover function.
        
        Returns
        -------
        skewness: float
            The value of the skewness for the scaling parameter used when executing the continuum_component_remover function."""

        # .flatten() makes a 2D numpy array into a 1D array.
        # .mean() when acted on a 2D numpy array would take the mean of each row or column.
        # Thus, to take the mean of all the entries in a 2D array, it must first be flattened.
        # .std() when acted on a 2D array would take the standard deviation of each row or column.
        # Thus, to take the standard deviation of all the entires in a 2D array, it must again first be flattened.
        # The result of "Residual_arr - Residual_arr.flatten().mean())/Residual_arr.flatten().std()" is a 2D numpy array.
        # .power(A, B) raises the power of each entry in A (a numpy array) to the power B (can either by a number of a numpy array)
        # .sum() takes the sum of all the entries in a row or column.
        # Hence, to take the sum of all the entries in the resultant 2D array, it needs to first be flattened.
        skewness = (1/(Residual_arr.size - 1)*(np.power((Residual_arr - Residual_arr.flatten().mean())/Residual_arr.flatten().std(), 3)).flatten().sum())

        return skewness


    def calcOptimalScaleFactor(self):
        self.skewness_vals = [] # Ensure skewness_vals is reset to an empty list

        for mu in self.factors:
            self.Residual_arr = Subtractor.getSubtractedImage(self.NB_arr, self.BB_arr, mu)
            skewness_value = Subtractor.skewnessCalculator(self.Residual_arr)
            self.skewness_vals.append(skewness_value)
            
        abs_skewness_vals = [] # Create an empty list to hold the magnitudes of the skewness values

        for i in range(len(self.skewness_vals)):
            abs_skewness_vals.append(abs(self.skewness_vals[i]))

        # The index of the skewness value which is closest to 0 is the same as the index of the minimum absolute skewness value
        index_optimal = abs_skewness_vals.index(min(abs_skewness_vals))

        optimal_factor = self.factors[index_optimal]

        return optimal_factor


    def plotPixelDistribution(self, scale_factor = 0.1, nbins = 301, lower_lim = -45000.0, upper_lim = 45000.0):
        """
        Plots a histogram of the pixel distribution for a specific scale factor.

        Parameters
        ----------
        scale_factor: float
            The scale factor to consider.
            Optional input.
            Default value is 0.1.
        
        nbins: integer
            Number of bins for the histogram.
            Optional input.
            Default value is 301.
        
        lower_lim: int, float
            Lower limit for the range of the x-axis for the histogram.
        
        upper_lim: int, float
            Upper limit for the range of the x-axis for the histogram.
        """

        self.Residual_arr = Subtractor.getSubtractedImage(self.NB_arr, self.BB_arr, scale_factor)

        # Flatten the array for histogram plotting
        flat_residual = self.Residual_arr.flatten()

        # Plot the histogram
        plt.hist(flat_residual, bins = nbins, range = (lower_lim, upper_lim))

        plt.text(0.5 * upper_lim, 0.9 * plt.ylim()[1], f"mu = {scale_factor}", ha='left', va='top')
        plt.text(0.5 * upper_lim, 0.4 * plt.ylim()[1], f"skew = {Subtractor.skewnessCalculator(self.Residual_arr): .2f}", ha='left', va='top')
        plt.yscale('log')

        # Label axes and display the plot
        plt.xlabel("Pixel Value")
        plt.ylabel("Frequency")
        plt.title("Pixel Distribution")
        plt.show()


    def plotSkewVsMu(self):
        plt.plot(self.factors, self.skewness_vals)

        # Set axis labels, anything in between $^\$ will be displayed as a mathematical symbol
        plt.xlabel('$^\mu$')
        plt.ylabel('s($^\mu$)')

    
    def getResultImage(self):
        # Turn Residual array back into a image class.
        pass