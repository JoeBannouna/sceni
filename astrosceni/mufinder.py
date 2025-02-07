import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astrosceni.image import Image
from scipy.interpolate import make_splrep, PPoly
from scipy.stats import skew

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
        print(self.mu_range)
        print(self.mu_resolution)
        print(((end - start)/self.mu_resolution))
        self.mu_linspace = np.linspace(start, end, int((end - start)/self.mu_resolution))
        self.is_paramaters_changed = True   # Flag that the parameters have changes

    def getSkewnessVals(self):
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

    def plotSkewnessVals(self):
        plt.xlabel('$^\mu$')
        plt.ylabel('s($^\mu$)')
        plt.plot(self.mu_linspace, self.getSkewnessVals())

    def getOptimalMus(self):
        skews = self.getSkewnessVals()
        print(np.sum(np.isnan(skews)))
        roots = PPoly.from_spline(make_splrep(self.mu_linspace, skews).derivative()).roots()
        roots = roots[roots < self.mu_range[1]]
        roots = roots[roots > self.mu_range[0]]
        return roots

    def getResultImages(self):
        optimal_mus = self.getOptimalMus()
        images = []
        for mu in optimal_mus:
            images.append(Image.subtract(self.narrow_band_image, self.broad_band_image, mu))
        return images