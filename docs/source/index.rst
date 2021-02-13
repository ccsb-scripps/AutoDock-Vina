AutoDock Vina: Molecular docking program
========================================

**AutoDock Vina** is one of the **fastest** and **most widely used** **open-source** docking engines. It is a turnkey computational docking program that is based on a simple scoring function and rapid gradient-optimization conformational search. It was originally designed and implemented by Dr. Oleg Trott in the Molecular Graphics Lab, and it is now being maintained and develop by the Forli Lab at The Scripps Research Institute.

Vina's design philosophy is not to require the user to understand its implementation details, tweak obscure search parameters, cluster results or know advanced algebra (quaternions). All that is required is the structures of the molecules being docked and the specification of the search space including the binding site. Calculating grid maps and assigning atom charges is not needed (when using Vina or Vinardo forcefields). Like in AutoDock 4, some receptor side chains can be chosen to be treated as flexible during docking. Additionally, it takes advantage of multiple CPUs or CPU cores on your system to significantly shorten its running time, which makes it faster than AutoDock 4 by *orders of magnitude*. Moreover, it significantly improves the average accuracy of the binding mode predictions compared to AutoDock 4.

Availability
------------
AutoDock Vina can be easily installed with all its dependencies using `pip` or `conda` package managers.

All **source code** is available under the `Apache License, version 2.0 <https://www.apache.org/licenses/LICENSE-2.0>`_ from `github.com/ccsb-scripps/AutoDock-Vina <https://github.com/ccsb-scripps/AutoDock-Vina>`_ and the Python Package index `pypi.org/project/Vina <pypi.org/project/Vina>`_.

Participating
-------------
Please report bugs or enhancement requests through the `Issue Tracker <https://github.com/ccsb-scripps/AutoDock-Vina/issues>`_.

AutoDock Vina is **open source** and welcomes your contributions. `Fork the repository on GitHub <https://github.com/ccsb-scripps/AutoDock-Vina>`_ and submit a pull request. Participate on the `developer mailing list <http://mgldev.scripps.edu/mailman/listinfo/autodock>`_.

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Manual

   introduction
   installation
   faq
   citations
   changes

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Tutorials

   docking_requirements
   docking_basic
   docking_flexible
   docking_multiple_ligands
   docking_zinc
   docking_hydrated
   docking_macrocycle
   docking_python

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Python Documentation

   vina

