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
