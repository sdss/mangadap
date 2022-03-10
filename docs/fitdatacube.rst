
.. include:: include/links.rst

.. _fitdatacube:

How to Fit a non-MaNGA Datacube
===============================

After the initial ingestion of the MaNGA datacube, the main
:ref:`execution-mangadap` is effectively independent of the specific data
format.  Although all algorithms were fine-tuned using MaNGA datacubes, the
pipeline can be executed on data from other instruments.  This document details
how this is done and gives an example using data from Keck's KCWI.

Overview
--------

The critical components that the DAP needs to perform its analysis modules are:

    #. spectra that are geometrically binned in wavelength,
    #. spectral errors, pixel masks, and an estimate of the spectral resolution,
    #. an estimate of the redshift that is accurate to within :math:`\pm 2000` km/s,
    #. conventions for how to identify the datacube being analyzed,
    #. a plan for the analyses to be performed, and
    #. a naming convention for the DAP output.

The first four components are handled by the datacube container class,
:class:`~mangadap.datacube.datacube.DataCube`.  Therefore, to run the DAP on a
custom datacube, you need to build a custom
:class:`~mangadap.datacube.datacube.DataCube` subclass.

The last two components are handled by the parameter class
:class:`~mangadap.config.analysisplan.AnalysisPlanSet`.  Depending on the
desired output path structure, it is likely that DAP analysis of a custom
datacube can use this class directly, instead of having to create a new
subclass.

.. _datacube_subclass:

Creating a new DataCube subclass
--------------------------------



  Instructions follow:

.. todo::
    - Primary attributes that must be defined are:
        - flux
        - ivar
        - mask
        - meta['z'] - Currently critical

    - Explain from_config



