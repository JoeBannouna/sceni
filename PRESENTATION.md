---
marp: true
theme: gaia
paginate: true
backgroundColor: #fff
---

![bg bottom:90% 90%](./images/original.png)
![bg bottom:90% 90%](./images/result.png)
# Subtraction of Continuum Emission from Narrow-band Images
<br><br><br><br>

<style scoped>small { font-size: 0.4rem; line-height: 0.1rem; }</style>

<small style="font-size: 10px;">Hong, Sungryong, et al. “Quantitative Method for the Optimal Subtraction of Continuum Emission from Narrow-Band Images: Skewness Transition Analysis.” Publications of the Astronomical Society of the Pacific, vol. 126, no. 935, 2014, pp. 79–99. JSTOR, https://doi.org/10.1086/674666. Accessed 10 Dec. 2024.</small>

---

<!-- _class: lead -->
# First Steps...

---

<!-- _class: lead -->
# Summarize paper

---

## Skewness
As $\mu$ increases, $NB - \mu \bullet BB$ transtion from undersubtracted to over subtracted
<br>

## skewness = $\frac{1}{N-1} \sum^N_{i=1} \left( \frac{x_i - m}{\sigma} \right)^3$

- Check animation

---

## One dimentional simulation

$I_{NB}(x) = E(x) + \sum S_i(x) + B_{NB} + Noise[\sigma_{NB}(x)]$

$I_{BB}(x) = \sum S_i'(x) + B_{BB} + Noise[\sigma_{BB}(x)]$

1. Extended line-emission $E(x)$
2. Emission from stars $\sum S_i(x)$
3. Background $B$
4. Combined noise of all kinds $σ(x)$

---


<!-- _class: lead -->
# Simulate NB & BB images and use POC to analyze results


---

## Proof of Concept
- One python/juptyer file.
- Show result image (NB image -  $\mu \bullet$ BB image) 
- Show result image for different $\mu$ scaling factors
- Each team member should try to code this example up

---

## Possible Ideas
<!-- Program should be able to: -->
- Superimpose stars from catalogue onto result image
- Mark all stars that have a mag larger than a set value
- Restrict the list of stars to within a given teh field of view
- Check if a given star vanished, check specific position
- Or list all stars that vanished (e.i. do not posses the emission line)
- Flexibility to choose catalogue (radio, etc..)
- Query a catalogue to find calibration sources (what is calibration source?)
- Show contours from BB image projected on (NB-μBB) image or vice versa?
<!-- - Calibration sources could be stars with periodicity < 2 -->
<!-- - hipparcos catalogue: for each star have ra, dec, periodicity, magnitude -->

---

<style scoped>section {font-size: 20px;}</style>

## What should program do? (in order)
- Find optimal $\mu$
  - User provides BB & NB images, and scale factor range (e.g. $\mu \in [0, 2]$)
  - User provides test region that will be used to calculate optimal $\mu$
  - Crop NB & BB data to the test region
  - Compute $NB-\mu \bullet BB$ for the range given in the test region
  - Calculate skewness of image for $\mu$ range
  - Find $\mu$ with skewness value closest to zero
- Remove stars
  - Get list of stars within the region of the full NB & BB images
  - Determine which stars can be seen (and thus need to be removed) (how?)
  - Remove the stars from NB & BB images (how?)
- Show final image
  - Subtract the edited NB & BB images using the optimal $\mu$
  - Show catalogue stars present in this region superimposed with user-filters
  - Have a feature to highlight catalogue stars that are not visible anymore (is this really necessary)
  - Display image (how? -> half pixels will be negative so need to think of how to plot) 


---

## Possible Classes
```python
class StarsFilter
  setCatalogue("name" or "http://link")
  getStarsInRegion(ra_1, ra_2, dec_1, dec_2, 
                  min_mag, max_mag, 
                  min_periodicity, max_periodicity)
```
```python
class StarsRemover
  setImage(image_header, imager_data)
  getStarsInImage()
  subtractStarsFromImage()
  getResultImage() # Decide if we need this later?
```
---

## Questions
- What is "Quality Criteria"
- In what way do we show contours in the images?
- How to/Should we calculate the fluxes?
- Selecting subset of difference image (NB-μBB) where stars dominate over nebulae for best accuracy, how??
  - <small>Make option for user to define this region/subset</small>
  - <small>Determine an objective method to find such subsets?</small>
    - <small><small>Use stars catalogue + nebulae catalogue to find subset with stars dominating over nebulae via some algorithm?</small></small>

---

## Project Management
Hello

---

## Codebase specifics
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
