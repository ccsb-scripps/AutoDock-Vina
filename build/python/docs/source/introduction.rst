AutoDock Vina
=============

Introduction
------------
AutoDock Vina (Vina) is one of the docking engines in the AutoDock Suite, together with AutoDock4, AutoDockGPU, AutoDockFR, and AutoDock-CrankPep, and arguably the most widely used docking engine. The reasons for this success are mostly due to its ease of use and its speed, compared to the other docking engines in the suite and elsewhere. Research groups around the world have built upon Vina-1.1.2 source code and improved the search algorithm (QuickVina),
made the interface more user friendly (Smina) and improved the scoring function [Vinardo].

Vina code is highly optimized and specialized, and one of its hallmarks is the very limited amount of user input necessary to perform a docking, as opposed to AutoDock4 (AD4), which requires a wide number of docking parameters to be defined for each calculation, including the explicit definition of the interaction grid maps to be used for each atom type. Although, performance and ease of use are achieved at the price of flexibility and adaptability. In fact, in AD4 the direct access to most of the engine internals can be exploited by simply tweaking the parameter files that are used to run the program (i.e., DPF, and force field parameter files) without source code modifications. Over the years, this has been successfully exploited to extend the software capabilities and add more functionalities, such as macrocyclic flexibility, metal coordination models, explicit waters, coarse grain models, and ligand irreversible bidning. 

Here we present the results of our exploratory effort in extending the Vina functionalities with some of the well-established methods developed for AD4, combining some of the functionalities and improvements of the two programs. 
