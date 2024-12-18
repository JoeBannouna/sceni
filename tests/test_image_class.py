def test_class_import():
  # if this throws an exception during loading, pytest will record a failure
  from astrosceni.image import Image

def test_image_load():
  from astrosceni.image import Image

  img1 = Image()
  img1.load("data/rim_Ha_wcs.fits")
  img1._printData()