# SCENI
**S**ubtraction of **C**ontinuum **E**mission from **N**arrow-band **I**mages.

## Description
- What does the package do?
- How does it do it? (breif high-level)
- Example images NB, BB, result image

## Installation
- Something like `pip install astrosceni` maybe

## Usage

Provide two images NB and BB
```python
import astrosceni

sceni = new astrosceni.sceni

sceni.NB('path/to/NB/image.fits')
sceni.BB('path/to/BB/image.fits')

sceni.setScaleFactorRange([0, 2])

# Find the scale factor
factor = sceni.findOptimalScaleFactor()

# find image difference
imageDiff = sceni.getImageDifference(scaleFactor=factor)

# automatically uses optimal factor if none specified
imageDiff = sceni.getImageDifference()

# Plot or save the image difference
imageDiff.plot()
imageDiff.save('path/to/save/image.fits')
```


```python
img = Image()
img.setNB('path/to/NB')
img.setBB('path/to/NB')

img.subtract()
img.addStarCircle(mag>2, 'blue')
img.addStarCircle(mag<=2, 'red')

img   # Show image (use __repr__)
```