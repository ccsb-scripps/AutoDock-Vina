Installation
============

Windows
-------

Vina is expected to work on Windows XP and newer systems.

**Installing**: 

Double-click the downloaded MSI file and follow the instructions

**Running**: 

Open the Command Prompt and, if you installed Vina in the default location, type

.. code-block:: bash

    "\Program Files\The Scripps Research Institute\Vina\vina.exe" --help

If you are using Cygwin, the above command would instead be

.. code-block:: bash

    /cygdrive/c/Program\ Files/The\ Scripps\ Research\ Institute/Vina/vina --help

See the Video Tutorial for details. Don't forget to check out Other Software for GUIs, etc.

Linux
-----
Vina is expected to work on x86 and compatible 64-bit Linux systems.

**Installing**: 

.. code-block:: bash

    tar xzvf autodock_vina_1_1_2_linux_x86.tgz

Optionally, you can copy the binary files where you want.

**Running**:

.. code-block:: bash

    ./autodock_vina_1_1_2_linux_x86/bin/vina --help

If the executable is in your PATH, you can just type "vina --help" instead. See the Video Tutorial for details. Don't forget to check out Other Software for GUIs, etc.

Mac
---

The 64 bit version is expected to work on Mac OS X 10.15 (Catalina) and newer. The 32 bit version of Vina is expected to work on Mac OS X from 10.4 (Tiger) through 10.14 (Mojave).

**Installing**:

.. code-block:: bash

    tar xzvf autodock_vina_1_1_2_mac_64bit.tgz   # 64 bit
    tar xzvf autodock_vina_1_1_2_mac.tgz         # 32 bit

Optionally, you can copy the binary files where you want.

**Running**:

.. code-block:: bash

    ./autodock_vina_1_1_2_mac_64bit/bin/vina --help     # 64 bit
    ./autodock_vina_1_1_2_mac/bin/vina --help           # 32 bit

If the executable is in your PATH, you can just type "vina --help" instead. See the Video Tutorial for details. Don't forget to check out Other Software for GUIs, etc.

Building from Source
--------------------

Attention: Building Vina from source is NOT meant to be done by regular users! (these instructions might be outdated)

- Step 1: **Install a C++ compiler suite**. On Windows, you may want to install Visual Studio; on OS X, Xcode; and on Linux, the GCC compiler suite.
- Step 2: **Install Boost**. Then, build and run one of the example programs, such as the Regex example, to confirm that you have completed this step. If you can't do this, please seek help from the Boost community.
- Step 3: **Build Vina**. If you are using Visual Studio, you may want to create three projects: lib, main and split, with the source code from the appropriate subdirectories. lib must be a library, that the other projects depend on, and main and split must be console applications. For optimal performance, remember to compile using the Release mode.

On OS X and Linux, you may want to navigate to the appropriate build subdirectory, customize the Makefile by setting the paths and the Boost version, and then type

.. code-block:: bash

    make depend
    make

