def test_class_import():
  # if this throws an exception during loading, pytest will record a failure
  from astrosceni.image import Image

def test_image_load():
  from astrosceni.image import Image

  img1 = Image()
  img1.load("data/rim_Ha_wcs.fits")
  img1._printData()

def test_images_crop():
  from astrosceni.image import Image

  img1 = Image()
  img1.load("data/rim_Ha_wcs.fits")
  img2 = Image()
  img2.load("data/rim_R_wcs.fits")

  img1._printInfo()
  img2._printInfo()

  img1.cropPixels(x_end=-20, y_start=20)
  img2.cropPixels(x_end=-20, y_start=20)

  img1._printInfo()
  img2._printInfo()
