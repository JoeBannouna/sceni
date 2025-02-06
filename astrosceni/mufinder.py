import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astrosceni.image import Image

class MuFinder:
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
            print(divisions)
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

    def getSkewnessVals(self):
        if (self.is_paramaters_changed):
            self.skewness_vals = []
            for mu in self.mu_linspace:
                data = (self.narrow_band_image.getImageData() - mu*self.broad_band_image.getImageData()).flatten()
                self.skewness_vals.append(np.sum(((data - np.mean(data)) / np.std(data, ddof=1)) ** 3) / (len(data) - 1))
        self.is_paramaters_changed = False
        return self.skewness_vals

    def plotSkewnessVals(self):
        plt.xlabel('$^\mu$')
        plt.ylabel('s($^\mu$)')
        plt.plot(self.factors, self.getSkewnessVals())

    def getOptimalMu(self):
        skews = self.getSkewnessVals()
        # BSPLINE CODE
        return 0.2

    def getResultImage(self):
        optimal_mu = self.getOptimalMu()
        return Image.subtract(self.narrow_band_image, self.broad_band_image, optimal_mu)