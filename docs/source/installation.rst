.. _installation:

Installation
============

.. warning::

    Currently, the python bindings and the binary installation are two separate processes. The installation of the python bindings does not include the Vina executable, and vice versa.

Unix- and Linux-based OS
------------------------

Vina is expected to work on x86 and compatible 64-bit Linux systems. The executable for the latest release are available here: `https://github.com/ccsb-scripps/AutoDock-Vina/releases <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_.

**Running**:

.. code-block:: bash

    ./vina_1.2.0_linux_x86_64 --help

If the executable is in your PATH, you can just type "vina --help" instead.

macOS
------

Vina is expected to work on macOS 10.15 (Catalina) and newer. The executable for the latest release are available here: `https://github.com/ccsb-scripps/AutoDock-Vina/releases <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_.

**Running**:

.. code-block:: bash

    ./vina_1.2.0_macos_x86_64 --help

If the executable is in your PATH, you can just type "vina --help" instead.

Python bindings
---------------

**AutoDock Vina installation using pip**:

.. note::

    When using ``pip``, it's good pratice to use a virtual environment and also the easiest solution (see ``meeko`` in :ref:`docking_requirements`). An example with the `Conda package manager <https://docs.conda.io/en/latest/>`_ is available further down.

.. code-block:: bash
    
    $ pip install -U numpy vina

**AutoDock Vina installation in a Conda environment**:

The Anaconda Python distribution, which can be download from `https://docs.continuum.io/anaconda/install <https://docs.continuum.io/anaconda/install/>`_. This is a Python distribution specially designed for scientific applications, with many of the most popular scientific packages preinstalled. Alternatively, you can use `Miniconda <https://conda.pydata.org/miniconda.html>`_, which includes only Python itself, plus the Conda package manager.

1. Begin by installing the most recent 64 bit, Python 3.x version of either Anaconda or Miniconda
2. Create a dedicated environment for ``AutoDock Vina``. This environment can be re-used for installing ``meeko`` (see :ref:`docking_requirements`):

.. code-block:: bash

    $ conda create -n vina python=3
    $ conda activate vina
    $ conda config --env --add channels conda-forge

3. And type the following command to install ``NumPy`` and ``AutoDock Vina``:

.. code-block:: bash

    $ conda install numpy
    $ pip install vina

Building from Source
--------------------

.. warning::

    Building Vina from source is NOT meant to be done by regular users!

- Step 1: **Install a C++ compiler suite**
    - Ubuntu/Debian: ``sudo apt-get install build-essentials``
    - macOS: Install Xcode from the `AppStore <https://apps.apple.com/fr/app/xcode/id497799835?mt=12>`_ and the Command Line Tools (CLT) from the terminal ``xcode-select --install``
- Step 2: **Install Boost and SWIG**
    - Ubuntu/Debian: ``sudo apt-get install libboost-all-dev swig``
    - macOS (with `Homebrew <https://brew.sh>`_): ``brew install boost swig``

- Step 3: **Build Vina**

    Start by downloading the lastest version of ``AutoDock Vina`` from github:

    .. code-block:: bash
    
        $ git clone https://github.com/ccsb-scripps/AutoDock-Vina

    To compile the binary (you might need to customize the Makefile by setting the paths to the Boost library):

    .. code-block:: bash

        $ cd AutoDock-Vina/build/linux/release
        $ make

    To compile the Python bindings:

    .. note::

        The ``Conda`` package manager is used here to easily install the several dependencies needed to build the ``Autodock-Vina`` python bindings (see above how to create a dedicated environment).

    .. code-block:: bash

        $ conda activate vina
        $ cd AutoDock-Vina/build/python
        $ conda install -c conda-forge numpy boost-cpp swig
        $ rm -rf build dist *.egg-info (to clean previous installation)
        $ python setup.py build install
