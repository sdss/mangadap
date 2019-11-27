Installation
============

Clone the repo
--------------

To download the DAP software and associated data, clone the `GitHub repo
<https://github.com/sdss/mangadap>`_ by executing:

    .. code-block:: bash

        git clone https://github.com/sdss/mangadap.git

This will create a ``mangadap`` directory in the current directory.
Although we try to keep the ``master`` branch of the repository stable,
we recommend using the most recent tag.  You can do so by executing:

    .. code-block:: bash

        cd mangadap
        ./checkout_current_tag


Install Python 3
----------------

The DAP is supported for Python 3 only.  To install Python, you can do
so along with a full package manager, like `Anaconda
<https://www.continuum.io/DOWNLOADS>`_, or you can install python 3
directly from `python.org <https://www.python.org/>`_.


Install the DAP code
--------------------

To install the DAP software, do one or more of the following (always
from within the top-level directory of the repo):

 * To perform an environment-level installation, run:

    .. code-block:: bash

        python3 setup.py install

 * On MacOSX, you may need to run:

    .. code-block:: bash

        CC=clang python3 setup.py install

 * To install the DAP in a way that changes you make to the repo are
   immediately available in your environment, run:

    .. code-block:: bash

        python3 setup.py develop

 * To install the DAP and ensure its dependencies are met, you can run:

    .. code-block:: bash

        pip3 install -e .

 * To install only the DAP dependencies, run:

    .. code-block:: bash

        pip3 install -r requirement.txt


Test your installation
----------------------

To test the installation, you can do one of the following:

 * Run the tests via the setup script:

    .. code-block:: bash

        python3 setup.py test

 * Run the tests using pytest directly:

    .. code-block:: bash

        cd python/mangdap/tests
        python3 -m pytest .


Local Environment Setup
-----------------------

The DAP uses environmental variable to define the paths to specific data
and other repositories.  If these are not defined, warnings will be
issued everytime the DAP is installed or imported.  The relevant
environmental variables, their default, and their usage are provided
below.

+----------------------------+-------------------------------------+------------------------------------------------+
|                   Variable |                             Default |                                       Comments |
+============================+=====================================+================================================+
| ``MANGADRP_VER``           | ``v2_4_3``                          | Version of the DRP, used for path construction |
+----------------------------+-------------------------------------+------------------------------------------------+
| ``MANGA_SPECTRO_REDUX``    | ``$HOME/MaNGA/redux``               | Root path for the reduced data                 |
+----------------------------+-------------------------------------+------------------------------------------------+
| ``MANGADAP_VER``           | ``mangadap.__version__``            | Version of the DAP, used for path construction |
+----------------------------+-------------------------------------+------------------------------------------------+
| ``MANGA_SPECTRO_ANALYSIS`` | ``$HOME/MaNGA/analysis``            | Root path for the analysis data                |
+----------------------------+-------------------------------------+------------------------------------------------+
| ``MANGACORE_VER``          | ``v1_6_2``                          | Version of MaNGA core (survey-level meta data) |
+----------------------------+-------------------------------------+------------------------------------------------+
| ``MANGACORE_DIR``          | ``$HOME/MaNGA/core/$MANGACORE_VER`` | Root path with the MaNGA core repository       |
+----------------------------+-------------------------------------+------------------------------------------------+

**Notes**

 * It is likely that you will never need to define ``$MANGACORE_VER``
   and ``$MANGACORE_DIR``.  These are only used in a specific mode of
   survey-level execution of the DAP.  See the ``rundap`` script.
 * The DAP expects to find the DRP ``LOGCUBE`` and ``LOGRSS`` files in
   the directory ``$MANGA_SPECTRO_REDUX/$MANGADRP_VER/[plate]/stack``,
   where ``[plate]`` is the desired plate number.  The ``LOGRSS`` files
   are required if you want to properly account for spatial covariance.
 * The DAP expects to find/write data to
   ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER``.
 * ``$MANGADAP_VER`` is only used to set the path names, not to select
   the specific version of the DAP that should be used.


These environmental variables can be added to, e.g., your
``.bash_profile`` file in your home directory or be included in a script
that is sourced when you want to run the DAP.  The added lines to your
`.bash_profile` file could look something like this:

    .. code-block:: bash

        export MANGA_SPECTRO_REDUX=/Volumes/MaNGA/redux
        export MANGA_SPECTRO_ANALYSIS=/Volumes/MaNGA/analysis

        export MANGADRP_VER=v2_4_3

        export MANGADAP_VER=2.2.1

        export MANGACORE_VER=v1_6_2
        export MANGACORE_DIR=$HOME/MaNGA/core/$MANGACORE_VER

**Note**: Some of these variables are also defined by `Marvin
<https://sdss-marvin.readthedocs.io/en/stable/installation.html>`_.  It
is possible to have both Marvin and the DAP point to the same directory,
but beware that this may mean that some of the files get overwritten!


Problems?
---------

We have limited support to offer installation help.  However, if you
have problems, particularly those that you think may be a more general
issue, please `submit an issue
<https://github.com/sdss/mangadap/issues>`_.

