from astropy.io import fits 
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import contour
import pandas as pd
from astroquery.vizier import Vizier

from astrosceni.image import Image

#Produces a contour overlay for a given image
class Contour:

    def __init__(self):
        self.contour_data = None
        self.contour_image = None

    #Generates a contour from a given image, and saves it
    def useImage(self, image):
        self.contour_data = np.array(image.current_image_data)

    #Plot the contour only
    def plot(self, custom_levels = 5):
        plt.figure()
        plt.contour(self.contour_data, levels = custom_levels, colors = 'black')
        plt.title("Contour Plot")
        plt.xlabel("x Axis")
        plt.ylabel("y Axis")
        plt.show()
        




        