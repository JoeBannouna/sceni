from astropy.io import fits 
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.visualization import simple_norm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import contour
import pandas as pd
from astroquery.vizier import Vizier
from scipy.ndimage import gaussian_filter

from astrosceni.image import Image

#Produces a contour overlay for a given image
class Contour:

    #class is initialized with no image
    def __init__(self):
        self.image_data = None

    #Load image into the program
    def useImage(self, image):
        self.image_data = image.getImageData()

    #Pre-processes data to enhance the visuals of faint object such as the hydrogen emission lines
    def preProcess(self):
        #
        smoothed_data = gaussian_filter(self.image_data, sigma = 2)
        # norm = simple_norm(smoothed_data, stretch='sqrt')
        print("Sigma = 3")
        return smoothed_data

    #Plot the contour only
    def plot(self):

        if np.all(self.image_data) == None:
            raise ValueError("No image loaded, load an image using useImage()")
        
        processed_data = self.preProcess()

        plt.figure(figsize=(10, 8))
        plt.contour(processed_data, levels = 10, cmap = 'viridis')
        plt.title("Contour plot")
        plt.xlabel("x-Axis")
        plt.ylabel("y-Axis")
        plt.colorbar(label = 'Intensity')
        plt.show()
        




        