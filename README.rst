ITKUltrasound
=============

.. image:: https://github.com/KitwareMedical/ITKUltrasound/workflows/Build,%20test,%20package/badge.svg

.. image:: https://img.shields.io/pypi/v/itk-ultrasound.svg
    :target: https://pypi.python.org/pypi/itk-ultrasound
    :alt: PyPI

.. image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
    :target: https://github.com/KitwareMedical/ITKUltrasound/blob/master/LICENSE)
    :alt: License

Purpose
=======

This module contains filters for use with the `Insight Toolkit (ITK)
<https://itk.org/>`_ that image formation and analysis of ultrasound
images.

Python Packages
===============

Python packages are available for Linux, macOS, and Windows. Install with::

  python -m pip install --upgrade pip
  python -m pip install itk-ultrasound

Build
=====

Build against ITK 5 or later. ITK should be configured with
*Module_SplitComponents*, *Module_Strain*, and *Module_BSplineGradient* set to *ON*.

References
==========

McCormick, M. An Open Source, Fast Ultrasound B-Mode Implementation for
Commodity Hardware. Insight Journal. 2010 January-June. URL:
https://hdl.handle.net/10380/3159

McCormick, M, Rubert, N and Varghese, T. Bayesian Regularization Applied to
Ultrasound Strain Imaging.  IEEE Transactions on Biomedical Engineering.
58 (6):1612-1620.  2011. PCMCID: PMC3092822.
https://doi.org/10.1109/TBME.2011.2106500
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3092822/

Aylward, S. R., McCormick, M. M., Kang H. J., Razzaque, S., R. Kwitt,
R., and M. Niethammer. Ultrasound spectroscopy. 2016 IEEE International
Symposium on Biomedical Imaging: From Nano to Macro, ISBI 2016 - Proceedings.
Prague, Czech Republic. 1013-1016. 2016.
https://dx.doi.org/10.1109/ISBI.2016.7493437
https://pdfs.semanticscholar.org/6bcd/1e7adbc24e15c928a7ad5af77bbd5da29c30.pdf
