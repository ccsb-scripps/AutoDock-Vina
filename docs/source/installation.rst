.. _installation:

Installation
============

Unix- and Linux-based OS
------------------------

Vina is expected to work on x86 and compatible 64-bit Linux systems.

**Installing**: 

.. code-block:: bash

    tar xzvf autodock_vina_1_2_0_linux_x86.tgz

Optionally, you can copy the binary files where you want.

**Running**:

.. code-block:: bash

    ./autodock_vina_1_2_0_linux_x86/bin/vina --help

If the executable is in your PATH, you can just type "vina --help" instead.

macOS
------

The 64 bit version is expected to work on macOS 10.15 (Catalina) and newer. The 32 bit version of Vina is expected to work on Mac OS X from 10.4 (Tiger) through 10.14 (Mojave).

**Installing**:

.. code-block:: bash

    tar xzvf autodock_vina_1_2_0_mac_64bit.tgz   # 64 bit
    tar xzvf autodock_vina_1_2_0_mac.tgz         # 32 bit

Optionally, you can copy the binary files where you want.

**Running**:

.. code-block:: bash

    ./autodock_vina_1_2_0_mac_64bit/bin/vina --help     # 64 bit
    ./autodock_vina_1_2_0_mac/bin/vina --help           # 32 bit

If the executable is in your PATH, you can just type "vina --help" instead.

Python bindings
---------------

**AutoDock Vina installation using Conda**:

The `Conda package manager <https://docs.conda.io/en/latest/>`_ is included as part of the Anaconda Python distribution, which can be download from `https://docs.continuum.io/anaconda/install <https://docs.continuum.io/anaconda/install/>`_. This is a Python distribution specially designed for scientific applications, with many of the most popular scientific packages preinstalled. Alternatively, you can use `Miniconda <https://conda.pydata.org/miniconda.html>`_, which includes only Python itself, plus the Conda package manager.

1. Begin by installing the most recent 64 bit, Python 3.x version of either Anaconda or Miniconda
2. (Optional, but highly suggested) If you want, you can create a dedicated environment for the `AutoDock Vina` package:

.. code-block:: bash

    $ conda create -n vina python=3
    $ conda activate vina
    $ conda install -c conda-forge numpy openbabel swig boost-cpp sphinx sphinx_rtd_theme
    $ conda install -c ccsb-scripps vina

3. And type the following command

.. code-block:: bash

    $ conda install -c ccsb-scripps vina

**AutoDock Vina installation using pip**:

.. code-block:: bash

    $ pip install vina

Building from Source
--------------------

.. warning::

    Building Vina from source is NOT meant to be done by regular users!

- Step 1: **Install a C++ compiler suite**
	- Ubuntu/Debian: ``sudo apt-get install build-essentials``
	- macOS: Install Xcode from the `AppStore <https://apps.apple.com/fr/app/xcode/id497799835?mt=12>`_ and the Command Line Tools (CLT) from the terminal ``xcode-select --install``
- Step 2: **Install Boost**
    - Ubuntu/Debian: ``sudo apt-get install libboost-all-dev``
    - macOS (with `Homebrew <https://brew.sh>`_): ``brew install boost``

- Step 3: **Build Vina**

    Start by downloading the lastest version of AutoDock Vina from github:

    .. code-block:: bash
    
        $ git clone https://github.com/ccsb-scripps/AutoDock-Vina

    To compile the binary (you might need to customize the Makefile by setting the paths to the Boost library):

    .. code-block:: bash

        $ cd AutoDock-Vina/build/linux/release
        $ make depend
        $ make

    To compile the Python bindings:

    .. code-block:: bash

        $ cd AutoDock-Vina/build/python
        $ python setup.py clean --all build install
