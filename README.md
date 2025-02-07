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
## Simplified Usage Example (more reccomended)

```python
from astrosceni import Image

img = Image(BB='path/to/BB', NB='path/to/NB')

# optionally
img.BB('path/to/BB')
img.NB('path/to/NB')

img.subtract()

img.showStars()
img.showContour()

img.addStars()
img.addContour()

img.setCatalogue()

img.addStarCircle(mag>2, color='blue')
img.addStarCircle(mag<=2, color='red')

img.addContour('file.fits', color=None)

img    # show image (use __repr__)
```

## TODO
- add scripts/ dir and inside a script to use the package and save the result image to a file from the command line
- use this package https://click.palletsprojects.com/en/stable/quickstart/#examples