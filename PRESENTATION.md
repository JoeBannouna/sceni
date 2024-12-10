---
marp: true
theme: gaia
paginate: true
backgroundColor: #fff
---

![bg bottom:90% 90%](./images/original.png)
![bg bottom:90% 90%](./images/result.png)
# Subtraction of Continuum Emission from Narrow-band Images

---

<!-- _class: lead -->
# First Steps...

---

## Proof of Concept
- One python/juptyer file.
- Show result image (NB image - mu $\bullet$ BB image) 
- Show result image for different (mu) scaling factors
- Each team member should try to code this example up

---

## Possible Ideas
Program should be able to:
- Superimpose stars from catalogue onto result image
- Mark all stars that have a mag larger than a set value
- Restrict the list of stars to within a given teh field of view
- Check if a given star vanished, check specific position
- Or list all stars that vanished (e.i. do not posses the emission line)
- Flexibility to choose catalogue (radio, etc..)
- Query a catalogue to find calibration sources (what is calibration source?)
<!-- - Calibration sources could be stars with periodicity < 2 -->
<!-- - hipparcos catalogue: for each star have ra, dec, periodicity, magnitude -->
- Show contours from BB image projected on (NB-BB) image or vice versa?

---

## Questions
- What is "Quality Criteria"
- In what way do we show contours in the images?
- How to/Should we calculate the fluxes?

---

## Project Management
Hello

---

## Codebase specifics
Testing