Installation
============

Clone the repo
--------------

To download the DAP software and associated data, clone the `GitHub repo
<https://github.com/sdss/mangadap>`_ by executing:

    .. code-block:: console

        git clone https://github.com/sdss/mangadap.git

This will create a ``mangadap`` directory in the current directory.
Although we try to keep the ``master`` branch of the repository stable,
we recommend using the most recent tag.  You can do so by executing:

    .. code-block:: console

        cd mangadap
        ./checkout_current_tag

.. warning::

    There is a distribution of the DAP that can be installed via
    `pip`, but the installation will be unsuccessful. We're still
    working out the bugs...

Install Python 3
----------------

The DAP is supported for Python 3 only.  To install Python, you can do
so along with a full package manager, like `Anaconda
<https://www.continuum.io/DOWNLOADS>`_, or you can install python 3
directly from `python.org <https://www.python.org/>`_.


Install the DAP code
--------------------

The preferred method to install the DAP and ensure its dependencies
are met is to, from the top-level, ``mangadap`` directory, run:

.. code-block:: console

    pip3 install -e .

This approach is preferred because it eases uninstalling the code:

.. code-block:: console
    
    pip3 uninstall sdss-mangadap

Alternatively, if you anticipate making changes to the DAP code, run:

.. code-block:: console

    python3 setup.py develop

----

To install only the DAP dependencies, run:

.. code-block:: console

    pip3 install -r requirements.txt


Test your installation
----------------------

To test the installation, you can do one of the following:

 * Run the tests via the setup script:

    .. code-block:: console

        python3 setup.py test

 * Run the tests using `pytest` directly:

    .. code-block:: console

        cd mangadap/tests
        python3 -m pytest .

Some tests requires a set of "remote" data that are not located in
the repo for space considerations. Downloading the data used by these
tests currently requires `SDSS Collaboration Access
<https://sdss-marvin.readthedocs.io/en/latest/installation.html#sdss-collaboration-access>`_.
The link in the last sentence points to a description of how this
access is granted for Marvin using a ``~\.netrc`` file. The DAP uses
the same ``~\.netrc`` file to authenticate access to the
``data.sdss.org`` host for downloading the test data. Once you have
your ``~\.netrc`` file, you can download the necessary test data and
rerun the tests to include usage of that data like this:

    .. code-block:: console

        python3 download_test_data.py
        cd mangadap/tests
        python3 -m pytest .


Local Environment Setup
-----------------------

The DAP uses environmental variables to define the paths to specific
data and other repositories. If these are not defined, warnings will
be issued every time the DAP is installed or imported. The relevant
environmental variables, their default, and their usage are provided
below.

+----------------------------+-------------------------------------+------------------------------------------------+
|                   Variable |                             Default |                                       Comments |
+============================+=====================================+================================================+
| ``MANGADRP_VER``           | ``v3_0_1`` (i.e., MPL-10)           | Version of the DRP, used for path construction |
+----------------------------+-------------------------------------+------------------------------------------------+
| ``MANGA_SPECTRO_REDUX``    | ``$HOME/MaNGA/redux``               | Root path for the reduced data                 |
+----------------------------+-------------------------------------+------------------------------------------------+
| ``MANGADAP_VER``           | ``mangadap.__version__``            | Version of the DAP, used for path construction |
+----------------------------+-------------------------------------+------------------------------------------------+
| ``MANGA_SPECTRO_ANALYSIS`` | ``$HOME/MaNGA/analysis``            | Root path for the analysis data                |
+----------------------------+-------------------------------------+------------------------------------------------+

These environmental variables can be added to, e.g., your
``.bash_profile`` file in your home directory or be included in a script
that is sourced when you want to run the DAP.  The lines added to your
``.bash_profile`` file could look something like this:

.. code-block:: bash

    export MANGA_SPECTRO_REDUX=/Volumes/MaNGA/redux
    export MANGADRP_VER=v3_0_1

    export MANGA_SPECTRO_ANALYSIS=/Volumes/MaNGA/analysis
    export MANGADAP_VER=3.0.1

.. note::

 * Importantly, note that ``$MANGADAP_VER`` is **only** used to set the
   path names, not to select the specific version of the DAP that
   should be used. The version of the DAP used is always the one
   installed by your python environment.
 * The DAP checks that these variables are defined *every time it is
   imported*. If they are not, warnings are raised and the defaults
   are used.
 * Some of these same variables are defined by `Marvin
   <https://sdss-marvin.readthedocs.io/en/stable/installation.html>`_.
   It is possible to have both Marvin and the DAP point to the same
   directory, but beware that this may mean that some of the files get
   overwritten!
 * Two additional variables (``$MANGACORE_VER`` and
   ``$MANGACORE_DIR``) are used in a specific mode of survey-level
   execution of the DAP. However, this is a niche usage mode and is
   effectively never used. See :ref:`execution-rundap`.
 * The DAP expects to find the DRP ``LOGCUBE`` *and* ``LOGRSS`` files
   in the directory
   ``$MANGA_SPECTRO_REDUX/$MANGADRP_VER/[PLATE]/stack``, where
   ``[PLATE]`` is the desired plate number. The ``LOGRSS`` files are
   required if you want to properly account for
   :ref:`spatialcovariance`. This path can be altered when executing
   the DAP.
 * The DAP expects to find/write data to
   ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER``. This path
   can be altered when executing the DAP, but the subdirectory
   structure used by the DAP to organize its outputs within this root
   directory cannot currently be changed.

Problems?
---------

We have limited support to offer installation help.  However, if you
have problems, particularly those that you think may be a more general
problem, please `submit an issue
<https://github.com/sdss/mangadap/issues>`_.

