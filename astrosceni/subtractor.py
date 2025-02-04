import numpy as np

import pandas as pd

import matplotlib.pyplot as plt

class Subtractor():

    def __init__(self, NB, BB):
        """
        Parameters:
        -----------
        NB: FITS file
            The narrowband FITS file.
            
        BB: FITS file
            The broadband FITS file.
        """
        self.NB_image = NB
        # convert NB FITS file to pandas dataframe
        self.NB_df = pd.DataFrame(self.NB_image.getImageData())

        self.BB_image = BB
        # convert BB FITS file to pandas dataframe
        self.BB_df = pd.DataFrame(self.BB_image.getImageData())

   
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

        if min >= max: raise ValueError('`min` must be smaller than `max`')
        if (max - min) < interval: raise ValueError('interval must be less than range')

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
            self.NB_image = self.NB_image[x_1: x_2, y_1: y_2]
            self.BB_image = self.BB_image[x_1: x_2, y_1: y_2]

    
    # define a function to do NB - mu*BB
    @staticmethod
    def continuum_component_remover(NB_df, BB_df, mu):
        """
        Calculates (Narrowband) - (scaling parameter)(Broadband).

        Parameters
        ----------
        NB_df: datafile 
            Narrowband datafile
        BB_df: datafile
            Broadband datafile
        
        mu: int, float
            The scale factor.
        
        Returns
        -------
        Residual: datafile
            The residual datafile"""

        Residual_df = NB_df - mu*BB_df

        return  Residual_df


    @staticmethod
    def skewness_calculator(Residual_df):
        """
        Calculates the skewness according to equation 1 from the paper (Sungryong Hong et al, 2014, Publications of the Astronomical Society of the Pacific).

        Parameters
        ----------
        Residual_df: dataframe
            The dataframe of the result from executing the continuum_component_remover function.
        
        Returns
        -------
        skewness: float
            The value of the skewness for the scaling parameter used when executing the continuum_component_remover function."""

        # N = Residual_df.size
        # Residual_df - Residual_df.mean())/Residual_df.std() calculates the difference between each value in the dataframe Residual_df and the mean of all values in the dataframe then divides each result by the standard deviation of the values in the dataframe.
        # .power(3) cubes each value from the above operation.
        # .values.sum() takes all the values from the resulting dataframe after the above operations and sums up all the values.
        skewness = (1/(Residual_df.size - 1)*((((Residual_df - Residual_df.values.mean())/Residual_df.values.std()).pow(3)).values.sum()))

        return skewness

    skewness_vals = []

    def calcOptimalScaleFactor(self):
        for mu in self.factors:
            self.Residual_df = Subtractor.continuum_component_remover(self.NB_df, self.BB_df, mu)
            skewness_vals = skewness_vals.append(Subtractor.skewness_calculator(self.Residual_df))
        for i in range(len(skewness_vals)):
            skewness_vals[i] = abs(skewness_vals[i])

        # The index of the skewness value which is closest to 0 is the same as the index of the minimum absolute skewness value
        index_optimal = skewness_vals.index(min(skewness_vals))

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

        self.Residual_df = Subtractor.continuum_component_remover(self.NB_df, self.BB_df, scale_factor)

        plt.hist(self.Residual_df, nbins, range = (lower_lim, upper_lim))
        plt.text(0, 1, f"mu = {scale_factor}", ha='left', va='top', transform=ax.transAxes)
        plt.text(0.2, 1, f"skew = {Subtractor.skewness_calculator(self.Residual_df)}", ha='left', va='top', transform=ax.transAxes)