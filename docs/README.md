# SCENI docs

## Table of contents
- [SCENI docs](#sceni-docs)
  - [Table of contents](#table-of-contents)
  - [Quick get started tutorial](#quick-get-started-tutorial)
  - [Classes](#classes)
    - [Introduction](#introduction)
    - [Image class](#image-class)
    - [StarsFilter class](#starsfilter-class)
    - [StarsRemover class](#starsremover-class)


## Quick get started tutorial

## Classes
Image

MuFinder

StarsFilter

### Introduction

### Image class
You can load a .fits file into an image class, where it will save the data and header info seperately. The image can be cropped through use of a cutout and is saved seperately as to not alter the originally saved data. Class is also used for plot and contour generation

### MuFinder class
Main purpose of this class is to determine what scaling factor (mu) is optimal in order to subtract the broadband image from the narrowband image and leave behind the emission spectrum being observed without losing too much information. Optimal mus are found through the calculation of how much each mu affects the final result (skewness).

### StarsFilter class
When given an image, uses an external star catalogue (by defauly Hipparcus) to identify stars within the bounds of the image. Gaussian fitting is also utilised to find the visible stars within the image. Finally these visible stars can be removed from the data once the image is passed in again. The simple usage class, filterStars(img), acts as all-encompassing function and provides all functionality with default arguments
