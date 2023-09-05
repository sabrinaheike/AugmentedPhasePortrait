
Version 2.0 (August 2023)

Made by Victoria Ralph.

Created under the supervision of Gail S.K. Wolkowicz and Sabrina H. Streipert with funding from NSERC (Undergraduate Student Research Award) through McMaster University.

Major Updates:
- Fixed bug in plotting the Y-nullcline which caused augmented_xy to take a long time to run.
- Changed cutoff function labels to correctly match the variable and value that points are mapped to.

- Plotting Next Iterative Operator Symbols 
    - Updated symbol plotting spacing so that symbols of next iterate operators would never overlap.
    - Removed symbols of next iterate operators at points mapped below the cutoff values (if given).
    - Added tolerance value so that symbols would not be plotted too close to root curves, preventing
      overlapping.

- Solving/Plotting Equilibria Points
    - Created equilibria.m function to numerically solve for nullcline intersections when it is not possible to        solve for them symbolically. After finding one point, the function recursively searches areas of the plane       excluding found points, so that all intersections are found.
    - Excluded equilibria points which were not real or where either of the X or Y functions were not defined.

 - Appearance of phase portrait
    - Updated colours of nullclines and cutoff functions (X-nullclines/root-curves/symbols of next iterate
      operators are plotted in shades of blue and Y-nullclines/root-curves/symbols of next iterate operators are       plotted in shades of red. Cutoff functions are plotted in shades of green.
    - Changed style of cutoff functions to dotted lines
    - Changed width of nullclines and varied nullcline widths so that overlapping nullclines would both be
      visible.



