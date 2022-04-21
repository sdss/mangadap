
.. include:: include/links.rst

Installation
============

Clone the repo
--------------

To download the DAP software and associated data, clone the `mangadap
GitHub repo`_ by executing:

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

The DAP is supported for Python 3 only. To install Python, you can do
so along with a full package manager, like `Anaconda`_, or you can
install python 3 directly from `python.org`_.


Install the DAP code
--------------------

The preferred method to install the DAP and ensure its dependencies
are met is to, from the top-level ``mangadap`` directory, run:

.. code-block:: console

    pip install -e .

This approach is preferred because it eases uninstalling the code:

.. code-block:: console
    
    pip uninstall sdss-mangadap

Alternatively, if you anticipate making changes to the DAP code, also install
the development dependencies:

.. code-block:: console

    pip install -e ".[dev]"

.. note::

    The use of the quotations above is shell dependent; e.g., you need them for
    zshell, but not for bash.  Also beware the exact characters used in the html
    above may not be the quotation characters you need for the command line
    (i.e., copy-pasting the line above may throw an error).

Test your installation
----------------------

To test the installation, make sure you have `pytest` installed and then

.. code-block:: console

    cd mangadap/tests
    pytest .

Some tests require a set of "remote" data that are not located in the repo for
space considerations. You can download these data by running the following
script:

.. code-block:: console

    python download_test_data.py

And then executing the tests with the same commands above.

Note that, if you have `SDSS Collaboration Access`_ and you have a `~/.netrc`
file (e.g., from using collaboration access with Marvin), this script will pull
the data from the (still public) `mangawork` directory on the SAS.


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
| ``MANGADRP_VER``           | ``v3_1_1`` (i.e., MPL-11)           | Version of the DRP, used for path construction |
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
    export MANGADRP_VER=v3_1_1

    export MANGA_SPECTRO_ANALYSIS=/Volumes/MaNGA/analysis
    export MANGADAP_VER=3.1.0

.. note::

 * Importantly, note that ``$MANGADAP_VER`` is **only** used to set the
   path names, not to select the specific version of the DAP that
   should be used. The version of the DAP used is always the one
   installed by your python environment.
 * The DAP checks that these variables are defined *every time it is
   imported*. If they are not, warnings are raised and the defaults
   are used.
 * Some of these same variables are defined by `Marvin`_. It is
   possible to have both Marvin and the DAP point to the same
   directory, but beware that this may mean that some of the files
   get overwritten!
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

If you have problems, particularly those that you think may be a more general
problem, please `Submit an issue`_.

