Standalone GUI: No need for MATLAB installed.

GUI uses the augmented_implicit.m and add_orbit.m version3.2 MATLAB functions for creating augmented phase portraits.

See Mac/Windows/Linux folders for dowloading instructions.

Version 3.2 Updates:
- Improvements to the nullcline combining algorithm to improve performance.
- New optional "Optimize Nullcline Combination Option".
  - This minimizes the total number of nullcline segments in the resulting phase portrait.
  - This is beneficial when there are multiple possible ways to combine nullcline segments.
  - It is not necessary to use the "Optimize Nullcline Combination Option" option in most cases.
    
User Guide:
Refer to the Readme file in the folder for each operating system (Mac, Windows, Linux) for operating system dependant downloading instructions.

User Guide:  

Using the app

-  Easily create augmented phase portraits using this desktop app.

- Input your X and Y equation after the "@(x,y)".

- Make sure equations are properly vectorized using a period (.) before *, /, and ^ operators. 

- Change the min and max values to adjust the axes limits.

- Change the accuracy using the slider. Increasing accuracy will increase the number of direction arrows and symbols of next iterate operators.

- The root curves and symbols of next iterate operators correspond only to the sections of nullclines calculated. When the next iterate operator is not able to be calculated a star is plotted in its place. To extend the root curves and symbols further, increase the nullcline calculation range by checking the box beside "Nullcline Calculation Range" and inputting your desired range of values. 

- To add thresholds indicated by green dash-dotted curves to show points (x,y) such that Xequation(x,y)=cutoffx and Yequation(x,y)=cutoffy check the box beside "Cutoff Values" and input the desired values into the boxes below. For example, to indicate where orbits leave the first quadrant, one would choose cutoffx = 0 and cutoffy = 0.

- The "Combine Nullcline Segments" option (suggested) joins nullcline segments which can form one continuous function. You can unclick this checkbox to deselect this feature.

- The "Optimize Nullcline Combination" option will join nullclines to mimimize the final number of nullclines for both X and Y.

- Add orbits by specifying the intial X and Y values and the number of iterations (the starting point counts as point 0). Then press the add orbit button.

- After adjusting all inputs, press the "Create Phase Portrait" button to plot your phase portrait.

Interpreting your phase portrait 

- The black arrows represent the direction field.

- X nullclines will be plotted with dashed curves in various shades of blue and Y nullclines in shades of red, and the associated root curves with solid curves of the same colour and shade as their nullcline.

- Symbols of next iterate operators will be in the same colour and shade as their corresponding nullcline. The symbols are triangles pointing left, right, up or down to represent the side of the nullcline the next iterate lies on. 

- In place of a triangle, a star will be plotted in the same colour and shade as the corresponding nullcline if the next iterate operator is unable to be evaluated. This is due to the point being mapped outside of the range of the nullcline. This occurs for nullclines that are function of x when the next iterate is not above/below the nullcline (the point is mapped too far to the left or right) so a value of up/down cannot be assigned to the next iterate operator. For nullclines that are functions of y this occurs when the next iterate is not to the left/right of the nullcline (the point is mapped too far up or down) so a value of left/right cannot be assigned to the next iterate operator.
