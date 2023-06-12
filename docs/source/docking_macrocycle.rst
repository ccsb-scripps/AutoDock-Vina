.. _macrocycle_docking:

Docking with macrocycles
========================

Macrocycles are made flexible by default since Meeko v0.3.0.
Thus, the basic docking tutorial suffices for docking macrocycles flexibly.

AutoDock Vina is not able to manage directly the flexibility associated with bonds in cyclic molecules. Different approaches can be used to dock macrocyclic molecules, like identifying one or more low energy conformations and docking them as different ligands, but generating them and docking them separately can be a time-consuming task. But AutoDock Vina (and AutoDock-GPU) has a specialized protocol to dock macrocycles while modeling their flexibility on-the-fly. 

The current implementation is described in our paper on the D3R Grand Challenge 4 (see below) and is summarized herein:

    1. One of the bonds in the ring structure is broken, resulting in an open form of the macrocycle that removes the need for correlated torsional variations, and enabling torsional degrees of freedom to be explored independently. 
    2. To each of the atoms of the broken bond, a dummy-atom is added. The distance between a dummy-atom and its parent atom is the same as the length of the broken bond, and the 1-3 angle also matches the original geometry. 
    3. During the docking, a linear attractive potential is applied on each dummy-atom to restore the bond resulting in the closed ring form. Thus, macrocycle conformations are sampled while adapting to the binding pocket, at the cost of increased search complexity with the added extra rotatable bonds. 

.. note::
    Even though dummy atoms are added and bonds are broken in the PDBQT representation of macrocycles,
    using Meeko's script mk_export.py (or class RDKitMolCreate from Python) will convert the docking
    output from PDBQT to MOL/SDF with correct connectivity and no dummy atoms. 


The followings papers can be cited if this method is important for your publication:

    - Holcomb, M., Santos-Martins, D., Tillack, A., & Forli, S. (2022). Performance evaluation of flexible macrocycle docking in AutoDock. QRB Discovery, 3, E18. doi:10.1017/qrd.2022.18
    - Forli, S., & Botta, M. (2007). Lennard-Jones potential and dummy atom settings to overcome the AUTODOCK limitation in treating flexible ring systems. Journal of chemical information and modeling, 47(4), 1481-1492
    - Santos-Martins, D., Eberhardt, J., Bianco, G., Solis-Vasquez, L., Ambrosio, F. A., Koch, A., & Forli, S. (2019). D3R Grand Challenge 4: prospective pose prediction of BACE1 ligands with AutoDock-GPU. Journal of Computer-Aided Molecular Design, 33(12), 1071-1081
