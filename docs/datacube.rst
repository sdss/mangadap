
.. _datacube:

Datacubes
=========

The majority of the DAP, in particular its full-spectrum fitting
modules, treat each spectrum independently (see :ref:`fitonespec`).
However, the primary input for the survey-level execution of the DAP
is a datacube.


The base class for all datacubes analyzed by the DAP is
:class:`~mangadap.datacube.datacube.DataCube`.  To run the MaNGA DAP on a
datacube from a currently unsupported instrument, you need to build a new 
:class:`~mangadap.datacube.datacube.DataCube` subclass.  Instructions follow:

.. todo::
    - Primary attributes that must be defined are:
        - flux
        - ivar
        - mask
        - meta['z'] - Currently critical

    - Explain from_config


