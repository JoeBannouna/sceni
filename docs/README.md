# SCENI docs
SCENI is a tool for processing astronomical images. Its main purpose is to load, alter and analyze '.fits' files, determine the optimal mu scaling factor for the subtraction of two images (NB and BB) to leave behind an emission spectrum and to remove visible stars utilising an externa catalogue

## Table of contents
- [SCENI docs](#sceni-docs)
  - [Table of contents](#table-of-contents)
  - [Quick get started tutorial](#quick-get-started-tutorial)
  - [Classes](#classes)
    - [Introduction](#introduction)
    - [Image class](#image-class)
    - [StarsFilter class](#starsfilter-class)
    - [StarsRemover class](#starsremover-class)
  - [Dependencies](#dependencies)

## Quick get started tutorial

## Classes
- Image
- MuFinder
- StarsFilter

### Introduction
Based on the following paper: https://www.jstor.org/stable/10.1086/674666

### Image class
The Image class handles the loading and manipulation of .fits files

**Key Methods**
- `load()`: Loads a '.fits' file.
- `cropPixels()`: Crops an image down to a specified amount of pixels
- `plot()`: Plots image with a large amount of user choice
- `subtract()`: Subtracts one image from another by a scaling factor mu

**Example:**
```python
nb_img = Image('example.fits')
nb_img.cropPixels(x_start, x_end, y_start, y_end)
nb_img.plot()
result_img = Image.subtract(nb_img, bb_imb, mu)
```

### MuFinder Class
Main purpose of this class is to determine what scaling factor (mu) is optimal in order to subtract the broadband image from the narrowband image and leave behind the emission spectrum being observed without losing too much information.

**Key Methods**
- `mufinder = MuFinder(nb_img, bb_img, mu_range = None)`: Loads both images into the MuFinder, if no mu range specified rough mu range will be generated
- `getOptimalMus()`: Returns optimal mu values in an array
- `getResultImages()`: Returns an array of pre-subtracted images with each optimal mu

**Example:**
```python
mufinder = MuFinder(nb_img, bb_img)
optimal_mus = mufinder.getOptimalMus()
result_imgs = mufinder.getResultImages()
```
### StarsFilter Class
When given an image, uses an external star catalogue (by defauly Hipparcus) to identify stars within the bounds of the image, determines which stars are visible and can subsequently remove them

**Key Methods**
- `setCatalogue()`: Sets catalogue to be hipparcus catalogue by default (can pass additional arguments to change this). Can also download catalogue to data folder if specified
- `setStarsInRegion(img)`: Returns dataframe of stars from the catalogue that are within the bounds of the iamge and their x and y pixel positions
- `setVisibleStars(img)`: Returns dataframe of stars from the catalogue that are within the bounds of the image only if they are considered visible by gaussian fitting
- `removeVisibleStars(img)`: Removes visible stars from image
- `filterStars(img)`: All encompassing function, includes all of the above

**Example:**
```python
starfilter.setCatalogue(download_catalogue = true)
starfilter.setStarsInRegion(img)
starfilter.setVisibleStars(img)
filtered_img = starfilter.removeVisibleStars(img)

#Alternatively, simpler use case
filtered_img = starfilter.filterStars(img)
```

## Dependencies
- Python 3.x
- NumPy
- SciPy
- AstroPy
- Matplotlib
- Pandas
- Astroquery
