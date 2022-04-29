
.. include:: include/links.rst

Installation
============

Install Python 3
----------------

The DAP is supported for Python 3 only. To install Python, you can do
so along with a full package manager, like `Anaconda`_, or you can
install python 3 directly from `python.org`_.

Python environment setup
------------------------

You are **strongly** encouraged to build an isolated python environment
specifically for the DAP.  This can be done using `virtualenv`_:

.. code-block:: console

    virtualenv dap
    source dap/bin/activate

or `conda`_:

.. code-block:: console

    conda create -n dap python=3.9
    conda activate dap

See the `Virtualenv documentation <https://virtualenv.pypa.io/en/latest/>`_
and/or `Managing Environments with Conda
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_
for more details. See also `virtualenvwrapper
<https://virtualenvwrapper.readthedocs.io/en/latest/>`_ as an option for more
easily managing `virtualenv`_ environments.

Install the DAP code
--------------------

You can install the DAP directly from `pip`_;

.. code-block:: console

    pip install sdss-mangadap

To uninstall:

.. code-block:: console
    
    pip uninstall sdss-mangadap

Test your installation
----------------------

To test the installation, make sure you have `pytest` installed and then

.. code-block:: console

    cd mangadap/tests
    pytest . -W ignore

Some tests require a set of "remote" data that are not located in the repo for
space considerations. You can download these data by running the following
script:

.. code-block:: console

    python download_test_data.py

And then executing the tests with the same commands above.

.. note::

    If you have `SDSS Collaboration Access`_ and you have a `~/.netrc` file
    (e.g., from using collaboration access with Marvin), this script will pull
    the data from the (still public) `mangawork` directory on the SAS.

Installation from source
------------------------

To install from the sources, first clone the `mangadap GitHub repo`_ by
executing:

.. code-block:: console

    git clone https://github.com/sdss/mangadap.git

This will create a ``mangadap`` directory in the current directory.
Although we try to keep the ``master`` branch of the repository stable,
we recommend using the most recent tag.  You can do so by executing:

.. code-block:: console

    cd mangadap
    ./checkout_current_tag

If you plan to develop the code (e.g., in your own fork of the code), you should
install the development dependencies as follows:

.. code-block:: console

    pip install -e ".[dev]"

.. note::

    The use of the quotations above is shell dependent; e.g., you need them for
    zshell, but not for bash.  Also beware the exact characters used in the html
    above may not be the quotation characters you need for the command line
    (i.e., copy-pasting the line above may throw an error).

Problems?
---------

If you have problems, particularly those that you think may be a more general
problem, please `Submit an issue`_.

