# Multiresolution Algorithms for Hybrid Imaging
This repository contains all code of the bachelor's thesis "Multiresolution Algorithms for Hybrid Imaging".
All of the Python code requires the scientific python suite, i.e. numpy, scipy, matplotlib.

As referenced in the thesis, we provide an animation of the energy density as an acoustical wave propagates through an illuminated optical medium:

<img src="https://i.imgur.com/vcWCTfL.gif" width="55%">

Also, we provide an animation (looping GIF, but also available as mp4 in this repository) of the multiscale image registration from the thesis:

<img src="https://i.imgur.com/D2gywgc.gif" width="55%">

The code is organized like this:
## `/multiscale_illustration`
`denoising.py` implements the TV minimization and denoising algorithm by Chambolle, as cited in the thesis. `multiscale_decomposition.py` may be used to compute a multiscale representation of a provided image, as demonstrated in the thesis. 
## `/fixed_point_algorithm`
C++ implementation of Ammari's fixed-point algorithm for ultrasound-modulated diffuse optical tomography. 
To compile this code, the finite-element framework deal.II is required. As mentioned in the thesis, the majority of this code 
was written in context of the seminar "Mathematics of Biomimetics" in the fall semester of 2020/2021.
## `/multiscale_image_registration`
Python implementation of the concrete multiscale image registration algorithm as described in the thesis. An mp4 file of the animation above is provided as well.
