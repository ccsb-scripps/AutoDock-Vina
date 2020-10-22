Docking with macrocycles
========================

Introduction
------------

AutoDock is not able to manage directly the flexibility associated with bonds in cyclic molecules, which leads to cyclic portions of the ligands to be considered as rigid. Different approaches can be used to dock macrocyclic molecules, like identifying one or more low energy conformations and docking them as different ligands, but generating them and docking them separately can be a time-consuming task. As an alternative, an indirect method may be used to manage the ring as a fully flexible entity and use the AutoDock conformation search to explore its flexibility. The method was initially developed for version 3.05, and now is implemented in version 4.2. The protocol converts the cyclic ligand into its corresponding acyclic form by removing a bond, and then docks the fully flexible molecule in the open form. A special atom type definition allows AutoDock to restore the original cycle structure during the calculation while exploring the cycle conformations during the search. The protocol can be subdivided in three main steps:

RING OPENING (a): by removing a bond, the ring is opened and the ligand is transformed to an acyclic form.

LIGAND PRE-PROCESSING (b):the ligand is processed following the standard AutoDockTools protocol, but the edge atoms are replaced with G atoms

DOCKING AND RING CLOSURE (c): the ligand is docked applying a linear attraction potential to the G-atoms that restore the cyclic structure

To restore the closed ring geometry a custom long range linear attraction potential is applied to these atoms during the docking calculation. This potential is effective at long range distances and guarantees the ring closure even with large cycles.

No extra maps are calculated for the G atoms because, for sake of evaluation of ligand-protein interaction, they are considered as normal carbon atoms. Therefore, C maps are used in their place. During the docking process, the potential guides the edge atoms next to each other resulting in an effective ring closure, while allowing the GA algorithm to explore the ring conformations.