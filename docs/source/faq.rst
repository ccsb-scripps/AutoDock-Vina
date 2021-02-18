Frequently Asked Questions
==========================

- **How accurate is AutoDock Vina?**

The predictive accuracy varies a lot depending on the target, so it makes sense to evaluate AutoDock Vina against your particular target first, if you have known actives, or a bound native ligand structure, before ordering compounds. While evaluating any docking engine in a retrospective virtual screen, it might make sense to select decoys of similar size, and perhaps other physical characteristics, to your known actives.

- **What is the difference between AutoDock Vina and AutoDock 4?**

AutoDock 4 (and previous versions) and AutoDock Vina were both developed in the Molecular Graphics Lab at The Scripps Research Institute. AutoDock Vina inherits some of the ideas and approaches of AutoDock 4, such as treating docking as a stochastic global opimization of the scoring function, precalculating grid maps (Vina does that internally), and some other implementation tricks, such as precalculating the interaction between every atom type pair at every distance. It also uses the same type of structure format (PDBQT) for maximum compatibility with auxiliary software.

However, the source code, the scoring funcion and the actual algorithms used are brand new, so it's more correct to think of AutoDock Vina as a new "generation" rather than "version" of AutoDock. The performance was compared in the original publication [*], and on average, AutoDock Vina did considerably better, both in speed and accuracy. However, for any given target, either program may provide a better result, even though AutoDock Vina is more likely to do so. This is due to the fact that the scoring functions are different, and both are inexact.

- **What is the difference between AutoDock Vina and AutoDock Tools?**

AutoDock Tools is a module within the MGL Tools software package specifically for generating input (PDBQT files) for AutoDock or Vina. It can also be used for viewing the results.

- **Can I dock two proteins with AutoDock Vina?**

You might be able to do that, but AutoDock Vina is designed only for receptor-ligand docking. There are better programs for protein-protein docking.

- **Will Vina run on my 64-bit machine?**

Yes. By design, modern 64-bit machines can run 32-bit binaries natively.

- **Why do I get "can not open conf.txt" error? The file exists!**

Oftentimes, file browsers hide the file extension, so while you think you have a file "conf.txt", it's actually called "conf.txt.txt". This setting can be changed in the control panel or system preferences.

You should also make sure that the file path you are providing is correct with respect to the directory (folder) you are in, e.g. if you are referring simply to conf.txt in the command line, make sure you are in the same directory (folder) as this file. You can use ls or dir commands on Linux/MacOS and Windows, respectively, to list the contents of your directory.

- **Why do I get "usage errors" when I try to follow the video tutorial?**

The command line options changed somewhat since the tutorial has been recorded. In particular, "--out" replaced "--all".

- **Vina runs well on my machine, but when I run it on my exotic Linux cluster, I get a "boost thread resource" error. Why?**

Your Linux cluster is [inadvertantly] configured in such a way as to disallow spawning threads. Therefore, Vina can not run. Contact your system administrator.

- **Why is my docked conformation different from what you get in the video tutorial?**

The docking algorithm is non-deterministic. Even though with this receptor-ligand pair, the minimum of the scoring function corresponds to the correct conformation, the docking algorithm sometimes fails to find it. Try several times and see for yourself. Note that the probability of failing to find the mininum may be different with a different system.

- **My docked conformation is the same, but my energies are different from what you get in the video tutorial. Why?**

The scoring function has changed since the tutorial was recorded, but only in the part that is independent of the conformation: the ligand-specific penalty for flexibility has changed.

- **Why do my results look weird in PyMOL?**

PDBQT is not a standard molecular structure format. The version of PyMOL used in the tutorial (0.99rc6) happens to display it well (because PDBQT is somewhat similar to PDB). This might not be the case for newer versions of PyMOL.

- **Any other way to view the results?**

You can also view PDBQT files in PMV (part of MGL Tools), or convert them into a different file format (e.g. using AutoDock Tools, or with "save as" in PMV)

- **How big should the search space be?**

As small as possible, but not smaller. The smaller the search space, the easier it is for the docking algorithm to explore it. On the other hand, it will not explore ligand and flexible side chain atom positions outside the search space. You should probably avoid search spaces bigger than 30 x 30 x 30 Angstrom, unless you also increase "--exhaustiveness".

- **Why am I seeing a warning about the search space volume being over 27000 Angstrom^3?**

This is probably because you intended to specify the search space sizes in "grid points" (0.375 Angstrom), as in AutoDock 4. The AutoDock Vina search space sizes are given in Angstroms instead. If you really intended to use an unusually large search space, you can ignore this warning, but note that the search algorithm's job may be harder. You may need to increase the value of the exhaustiveness to make up for it. This will lead to longer run time.

- **The bound conformation looks reasonable, except for the hydrogens. Why?**

AutoDock Vina actually uses a united-atom scoring function, i.e. one that involves only the heavy atoms. Therefore, the positions of the hydrogens in the output are arbitrary. The hydrogens in the input file are used to decide which atoms can be hydrogen bond donors or acceptors though, so the correct protonation of the input structures is still important.

- **What does "exhaustiveness" really control, under the hood?**

In the current implementation, the docking calculation consists of a number of independent runs, starting from random conformations. Each of these runs consists of a number of sequential steps. Each step involves a random perturbation of the conformation followed by a local optimization (using the Broyden-Fletcher-Goldfarb-Shanno algorithm) and a selection in which the step is either accepted or not. Each local optimization involves many evaluations of the scoring function as well as its derivatives in the position-orientation-torsions coordinates. The number of evaluations in a local optimization is guided by convergence and other criteria. The number of steps in a run is determined heuristically, depending on the size and flexibility of the ligand and the flexible side chains. However, the number of runs is set by the exhaustiveness parameter. Since the individual runs are executed in parallel, where appropriate, exhaustiveness also limits the parallelism. Unlike in AutoDock 4, in AutoDock Vina, each run can produce several results: promising intermediate results are remembered. These are merged, refined, clustered and sorted automatically to produce the final result.

- **Why do I not get the correct bound conformation?**

It can be any of a number of things:

    1. If you are coming from AutoDock 4, a very common mistake is to specify the search space in "points" (0.375 Angstrom), instead of Angstroms.
    2. Your ligand or receptor might not have been correctly protonated.
    3. Bad luck (the search algorithm could have found the correct conformation with good probability, but was simply unlucky). Try again with a different seed.
    4. The minimum of the scoring function correponds to the correct conformation, but the search algorithm has trouble finding it. In this case, higher exhaustiveness or smaller search space should help.
    5. The minimum of the scoring function simply is not where the correct conformation is. Trying over and over again will not help, but may occasionally give the right answer if two wrongs (inexact search and scoring) make a right. Docking is an approximate approach.
    6. Related to the above, the culprit may also be the quality of the X-ray or NMR receptor structure.
    7. If you are not doing redocking, i.e. using the correct induced fit shape of the receptor, perhaps the induced fit effects are large enough to affect the outcome of the docking experiment.
    8. The rings can only be rigid during docking. Perhaps they have the wrong conformation, affecting the outcome.
    9. You are using a 2D (flat) ligand as input.
    10. The actual bound conformation of the ligand may occasionally be different from what the X-ray or NMR structure shows.
    11. Other problems 

- **How can I tweak the scoring function?**

You can change the weights easily, by specifying them in the configuration file, or in the command line. For example

    ``vina --weight_hydrogen -1.2 ...``

doubles the strenth of all hydrogen bonds.

- **Functionality that would allow the users to create new atom and pseudo-atom types, and specify their own interaction functions is planned for the future.**

This should make it easier to adapt the scoring function to specific targets, model covalent docking and macro-cycle flexibility, experiment with new scoring functions, and, using pseudo-atoms, create directional interaction models.

Stay tuned to the AutoDock mailing list, if you wish to be notified of any beta-test releases.

- **Why don't I get as many binding modes as I specify with "--num_modes"?**

This option specifies the maximum number of binding modes to output. The docking algorithm may find fewer "interesting" binding modes internally. The number of binding modes in the output is also limited by the "energy_range", which you may want to increase.

- **Why don't the results change when I change the partial charges?**

AutoDock Vina ignores the user-supplied partial charges. It has its own way of dealing with the electrostatic interactions through the hydrophobic and the hydrogen bonding terms. See the original publication [*] for details of the scoring function.

- **I changed something, and now the docking results are different. Why?**

Firstly, had you not changed anything, some results could have been different anyway, due to the non-deterministic nature of the search algorithm. Exact reproducibility can be assured by supplying the same random seed to both calculations, but only if all other inputs and parameters are the same as well. Even minor changes to the input can have an effect similar to a new random seed. What does make sense discussing are the statistical properties of the calculations: e.g. "with the new protonation state, Vina is much less likely to find the correct docked conformation".

- **How do I use flexible side chains?**

You split the receptor into two parts: rigid and flexible, with the latter represented somewhat similarly to how the ligand is represented. See the section "Flexible Receptor PDBQT Files" of the AutoDock4.2 User Guide (page 14) for how to do this in AutoDock Tools. Then, you can issue this command: vina --config conf --receptor rigid.pdbqt --flex side_chains.pdbqt --ligand ligand.pdbqt. Also see this write-up on this subject.

- **How do I do virtual screening?**

Please see the relevant section of the manual. Please note that a variety of docking management applications exist to assist you in this task.

- **I have ideas for new features and other suggestions.**

For proposed new features, we like there to be a wide consensus, resulting from a public discussion, regarding their necessity. Please consider starting or joining a discussion on the AutoDock mailing list.

- **Will you answer my questions about Vina if I email or call you?**

No. Vina is community-supported. There is no obligation on the authors to help others with their projects. Please see this page for how to get help. 