from astrosceni.image import Image
from astrosceni.mufinder import MuFinder
from astrosceni.stars_filter import StarsFilter

# This file serves as an example for simple usage of the library

# Load data into two image objects
nb_img = Image("data/bs_h_ave_wcs.fits")
bb_img = Image("data/bs_r_ave_wcs.fits")

# Crop images 100 pixels in from all sides to remove camera artifacts around edges
nb_img.cropPixels(100, -100, 100, -100)
bb_img.cropPixels(100, -100, 100, -100)

# Create a MuFinder object, load in cropped images
mufinder = MuFinder(nb_img, bb_img, mu_resolution = 0.05)

# Obtain optimal mu value from mufinder
optimal_mus = mufinder.getOptimalMus()

# Initialize a starfilter object and remove visible stars from images (stars by default are found through use of hipparcus catalogue)
starsFilter = StarsFilter()
filtered_nb = starsFilter.filterStars(nb_img)
filtered_bb = starsFilter.filterStars(bb_img)

# Subtract bb from nb by a scaling factor of the optimal mu
result_img = Image.subtract(filtered_nb, filtered_bb, optimal_mus[0])

# Plot image
result_img.plot(cmap = 'viridis') 