# Fluorescence Microscope Utilities
During my Postdoc with [Pasteur Institute](https://pasteur.fr) in Paris, I contributed to the open source project called [Icy](https://icy.bioimageanalysis.org/) for biological image analysis. The lab is lead by [Dr. Jean-Christophe Olivo-Marin](https://research.pasteur.fr/en/member/jean-christophe-olivo-marin/)

This repository contains 3 tools:
1. [Nyquist Sampling](https://github.com/praveenpankaj/Fluorescence-Microscope-PSF/blob/master/NyquistSampling.jar) 
This plugin calculates the Nyquist sampling, in the radial and axial direction for your Microscope. 

2. [Widefield PSF](https://github.com/praveenpankaj/Fluorescence-Microscope-PSF/blob/master/WideFieldPSF.jar)
This plugin generates the point-spread function (PSF) for a Widefield Fluorescence Microscope using a scalar diffraction-limited model (Stokseth). It can generate PSF with spherical aberration if there is refractive index mismatch between objective immersion and specimen mounting media.

3. [Widefield Macroscope PSF](https://github.com/praveenpankaj/Fluorescence-Microscope-PSF/blob/master/MacroscopePSF.jar)
This plugin generates the point-spread function (PSF) for a Widefield Fluorescence Macroscope. Due to the telescopic system, the PSF is not spatially invariant.


## References: 
[1] Sheppard, C.J.R. The spatial frequency cut-off in three-dimensional imaging. (1986a). Optik 72 No. 4 131-133. 
[2] Sheppard, C.J.R.The spatial frequency cut-off in three-dimensional imaging II. (1986b). Optik 74 No. 3, pp. 128-129. 
[3] https://svi.nl/tiki-index.php?page=NyquistRate
[4] P. A. Stokseth (1969). `Properties of a defocused optical system'. J. Opt. Soc. Am. A 59:1314-1321. 
[5] P. Pankajakshan (2009). Blind Deconvolution for Confocal Laser Scanning Microscopy. Ph.D. thesis, Universite de Nice Sophia-Antipolis.
[6] P. Pankajakshan et al., "Point-spread function model for fluorescence MACROscopy imaging," 2010 Conference Record of the Forty Fourth Asilomar Conference on Signals, Systems and Computers, 2010, pp. 1364-1368, doi: 10.1109/ACSSC.2010.5757756.
