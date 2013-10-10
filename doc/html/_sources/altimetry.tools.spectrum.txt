.. _spectral_tools:

Spectral analysis
=================

The :mod:`altimetry.tools.spectrum` module contains tools dedicated to spectral analysis.

About spectral analysis
-----------------------

Spectral analysis of along-track data is a common thing. There are 2 main steps when computing a spectrum:
  
  * **preprocess the data**
     
    It consists in detecting gaps, interpolating over short gaps and rejecting longer gaps, subsampling the data into subsegments of valid data of a given length.
    
    This step is performed using :func:`altimetry.tools.spectrum.preprocess`
     
  * **compute the spectrum**
  	
    This step is made through a transform of the signal to the spectral domain (eg. FFT). Then frequency, energy and power spectral densities are computed and averaged. It is also possible to use **spectral tapers** to lower the noise of the spectrum. 

	This step is performed using :func:`altimetry.tools.spectrum.spectral_analysis` (and :func:`altimetry.tools.spectrum.get_spec` at lower level)

Notes on spectral tapering
++++++++++++++++++++++++++

Tapering and padding are mathematical manipulations sometimes performed on the time series before periodogram analysis to improve the statistical properties of the spectral estimates or to speed up the computations.

**Tapering can be applied:**
  * to reduce the noise level by oversampling the data in overlapping subsegments (eg. when we don't have enough samples)
  * to better localise spectral peaks and changes in the spectral slope.

**However, you should be aware that**:
  * tapering may induce a loss of overall energy, resulting the tapered spectrum to be under (though less noisy) the original spectrum.
  * oversampling the data will result in removing a part of the lower frequencies because of the shorter subsegments.

:func:`altimetry.tools.spectrum.preprocess` allows using tapers through its ``tapering`` keyword.

.. warning:: though it is taken into account in :func:`altimetry.tools.spectrum.spectral_analysis`, energy loss caused by the tapering may not be properly resolved.
  
  It may be therefore necessary to correct this loss by multiplying the tapered spectrum by the ratio of energies of both spectra :math:`\frac{E_{original}}{E_{tapered}}` 

Notes on AR spectrum (auto-regression methods)
++++++++++++++++++++++++++++++++++++++++++++++

AR (auto-regressive methods) can be used to model a spectrum from the signal.

Such method, as the **Yule-Walker equations**, can be used to model the spectrum, and therefore:
  * clean the spectrum (by having an auto-regression approach)
  * compute the energy (or power) at any frequency (ie. not being dependant on the length of input array).

This approach is made possible through the :keyword:`ARspec` keyword of :func:`altimetry.tools.spectrum.spectral_analysis` (itself calling :func:`altimetry.tools.spectrum.yule_walker_regression`).

List of useful functions
------------------------

  * :func:`altimetry.tools.spectrum.spectral_analysis` : Compute the average spectrum over a set of data
  * :func:`altimetry.tools.spectrum.preprocess` : Preprocess the data to be admissible to spectral analysis
  * :func:`altimetry.tools.spectrum.get_slope` : Compute the spectral slope
  * :func:`altimetry.tools.spectrum.optimal_AR_spectrum` : Get the order of the optimal AR spectrum 

Functions
---------

.. automodule:: altimetry.tools.spectrum
	:members: spectral_analysis, preprocess, get_kx, get_spec, get_segment, get_slope, yule_walker, yule_walker_regression, optimal_AR_spectrum