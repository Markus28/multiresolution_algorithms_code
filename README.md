# Multiresolution Algorithms for Hybrid Imaging
This repository contains all code of the bachelor's thesis "Multiresolution Algorithms for Hybrid Imaging".
All of the Python code requires the scientific python suite, i.e. numpy, scipy, matplotlib.

As referenced in the thesis, we provide an animation of the energy intensity as an acoustical wave propagates through an illuminated optical medium:

<img src="https://i.imgur.com/vcWCTfL.gif" width="50%">

Also, we provide an animation of the multiscale image registration from the thesis:

<img src="https://i.imgur.com/N0KGhFN.gif" width="50%">

The code is organized like this:
## `/multiscale_illustration`
Python implementation of Chambolle's TV minimization algorithm, as cited in the thesis. `denoising.py` implements the TV minimization and denoising algorithm by Chambolle. 
`multiscale_decomposition.py` can be used to compute a multiscale representation of a provided image, as demonstrated in the thesis. 
## `/fixed_point_algorithm`
C++ implementation of Ammari's fixed-point algorithm for ultrasound-modulated diffuse optical tomography. 
To compile this code, the finite-element framework deal.II is required. As mentioned in the thesis, the majority of this code 
was written in context of the seminar "Mathematics of Biomimetics" in the fall semester of 2020.
## `/multiscale_image_registration`
Python implementation of the concrete multiscale image registration algorithm as described in the thesis. An mp4 file of the animation above is provided as well.
