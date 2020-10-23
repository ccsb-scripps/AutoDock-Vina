Changes
=======

1.2.0
-----
  
  - Refactoring of the Vina code to be able to be used as a library
  - Can load external maps
  - Support of simultaneous docking of multiple ligands and batch mode for virtual screening
  - Support of macrocycle molecules
  - Addition of the hydrated docking protocol
  - Added AutoDock4 and Vinardo forcefields
  - Added Python bindings for Python 3

1.1.2
-----
  
  - Bug fix: the affinities reported with flexible side chains were
    sometimes incorrect. (This problem had no effect on the predicted
    binding modes or their relative ranking, or on rigid receptor
    docking)

  - Bug fix: the individual term values, before weighting, reported
    only with "score_only", were often wrong. (This problem
    had no effect on the predicted binding modes or affinities)

  - The predicted binding modes are made less redundant by requiring
    that no two reported modes differ by less than 1 Angstrom RMSD(U.B.),
    counting all movable heavy atoms, including those in the side chains.
    Previously, Vina avoided this kind of redundancy during the actual
    docking, but made no such guarantee w.r.t. the output because
    of the subsequent refinement stage that could move different binding
    modes closer

  - Bug fix: in some very unusual cases, numerical rounding errors would 
    accumulate, leading, internally, to a distorted ligand structure.
    The self-checks in Vina would catch this and make the program quit
    rather than write invalid output. This numerical error accumulation 
    is believed to have been fixed in this update

  - A warning is added when the search space is greater than
    27000 Angstrom^3 in volume, as confusing Angstroms with
    "grid points" is a frequent user error

1.1.1
-----

  - Source code released

  - License changed to Apache

  - Extensive code clean-up that however should not have any effect
    on the end users

  - Unix binary packages changed to simple archives

1.1.0
-----

  - Added "advanced options" intended to be primarily used by people 
    interested in methods development.

    These should allow scoring without minimization; performing 
    local optimization only; randomizing the input with no search 
    (this is useful for testing docking software); changing the 
    weights from their default values; displaying the individual 
    contributions to the intermolecular score, before weighting 
    (these are shown with "--score_only") 
  
1.0.3
-----

  - Mac OS X 10.6 (Snow Leopard) support added

  - A bug related to (very rarely) not finding any modes when the
    "num_modes" setting is low (e.g. 1) has been fixed.

1.0.2
-----

  - Fixed a bug that excluded conformations in which only a hydrogen atom was
    outside the search space.

1.0.1
-----

  - When the search space is too small or misplaced, so that the ligand and 
    flexible side chain conformations can not be found within it, Vina will
    no longer produce any "binding modes" and print a warning instead

  - Minor improvements in user-friendliness of the error messages.

1.0
---

  - Support for earlier versions of Linux added, such as Debian 4.0

  - Friendlier and more specific error messages for PDBQT syntax errors

  - "all" parameter is now called "out" and is always enabled. Individual 
    binding modes are no longer written separately. If needed, the multimodel
    output can be split into individual models using a separate program
    called "vina_split", included in the distribution

  - The default "num_modes" value changed to 9.
