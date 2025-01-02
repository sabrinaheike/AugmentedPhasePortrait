Desktop app and Matlab code for creating the augmented phase portrait for discrete planar maps and adding orbits. The GUIs are standalone in 3 flavors, Mac, Windows, and Linux, and they do not require Matlab to be installed.
However, if you do have Matlab installed and prefer not to download .exe files, you can create the GUI by downloading the ascii files: 2 each with extensions .m and .png, and 1 each with extensons .pjr and .mlapp. Then click on the .pjr file and follow the instructions and the AugmentedPhasePlane.exe file will be created. It will prompt you to download the free Matlab Runtime program that is compatible with your version of Matalab.

Improvements from version3.1 code:
- Fixed bug in assigning nullclines as functions of x and y to ensure that next iterate operator values and root curve points are calculated properly.
- Improved nullcline combingin algorithm. Increased minimum distance between nullcline segment endpoints considered close enough to combine segments into a single nullcline.

Augmented phase portraits introduced in S.H. Streipert, G.S.K. Wolkowicz: An augmented phase plane approach for discrete planar maps: Introducing next-iterate operators, https://doi.org/10.1016/j.mbs.2022.108924

Created by: Victoria Ralph, ralphv1@mcmaster.ca, McMaster University.
Supervised by Gail S. K. Wolkowicz, McMaster University and Sabrina H. Streipert, University of Pittsburgh.
