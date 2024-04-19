MATLAB program to create augmented phase portraits.

Improvements from version3.0 code: 
-  Add orbits to your phase portrait using the add_orbit.m function
-  New option to automatically combine nullcline segments using the CombineNullcline argument
-  Removal of erroneous nullclines which appear due to discontinuities in the function

Instructions for using augmented_implicit.m: 
1. Download file and open in MATLAB
2. Call function using augmented_implicit(inputs)

Required Inputs: (f,g,minx, maxx, miny, maxy)

        f = @(x,y) f(x,y); right-hand side of the X-equation f = f(x,y)
        g = @(x,y) g(x,y); right-hand side of the Y-equation g = g(x,y) 
      (Make sure f and g equations are properly vectorized using a period (.) before *, /, and ^ operators)
        minx = number; minimum x-value to be considered for plotting
        maxx = number; maximum x-value to be considered for plotting
        miny = number; minimum y-value to be considered for plotting
        maxy = number; maximum y-value to be considered for plotting

Optional Inputs: add as Name-Value pairs
    Accuracy = number; optional parameter for number of arrows and signs of 
       next-iterate operators; default value = 20
    CutoffValues = 1x2 double; if numbers given then dash-dotted curves are
       added in plot to show points (x,y) such that f(x,y)=cutoffx and
       g(x,y)=cutoffy
    PortraitAxes = axes object; specifies the axes to plot the augmented
       phase portrait
    NullclineRange = 1x4 double; Optional range over which to calculate the
       nullclines. The root curves and symbols of next iterate operators
       correspond only to the sections of nullclines that are plotted. In order
       to make the root curves and symbols extend further, add a nullcline calculation 
       range which is larger than the minimum and maximum values for plotting.
       The nullcline calculation values are to be inputed in the format
       [minx, maxx, miny, maxy].
    CombineNullcline = string "on" or "off"; The automatic setting is for the CombineNullcline option is "on". When on, the program will     
      automatically sort and attempt to combine nullcline segments which can be joined to make a single continuous function. This reduces the       number of nullclines and symbols in the portrait. 

Interpreting your phase portrait 

- The black arrows represent the direction field.

- X nullclines will be plotted with dashed curves in various shades of blue and Y nullclines in shades of red, 
    and the associated root curves with solid curves of the same colour and shade as their nullcline.

- Symbols of next iterate operators will be in the same colour and shade as their corresponding nullcline. 
    The symbols are triangles pointing left, right, up or down to represent the side of the nullcline the next iterate lies on. 

- In place of a triangle, a star will be plotted in the same colour and shade as the corresponding nullcline if the next iterate operator is unable to be evaluated. This is due to the point being mapped outside of the range of the nullcline. This occurs for nullclines that are function of x when the next iterate is not above/below the nullcline (the point is mapped too far to the left or right) so a value of up/down cannot be assigned to the next iterate operator. For nullclines that are functions of y this occurs when the next iterate is not to the left/right of the nullcline (the point is mapped too far up or down) so a value of left/right cannot be assigned to the next iterate operator.

Instructions for using add_orbit.m: 
  1. Download file and open in Matlab
  2. Add orbits to your phase portait using the command add_orbit(args)

Required Inputs: (f,g,Xo,Yo,iterations)
  f = @(x,y) f(x,y); right-hand side of the X-equation f = f(x,y (the same as for the phase portrait)  
  g = @(x,y) g(x,y); right-hand side of the Y-equation g = g(x,y) (the same as for the phase portrait)  
  Xo = x-coordinate of initial conditions for orbit  
  Yo = y-coordinate of initial conditions for orbit  
  iterations = The number of iterations for your orbit (the initial point is point 0)  

Optional Input: add as Name-Value pair
  OrbitAxes = axes object; specify the axes on which to plot the orbit


Example: 
To get the augmented phase portrait for the system: 
x_(t+1) = 2x_t/(1+x_t+0.3y_t),
y_(t+1) = 3y_t/(1+2y_t+0.6x_t)

plotted in [0,2]x[0,3] with
accuracy 30,
cutoff values of 0 for both x and y, and 
using a nullcline calculation range of [0,4]x[0,5]

Make sure equations are properly vectorized using a period (.) before *, /, and ^ operators

Type:

> augmented_implicit(2 * x ./ (1 + x + 0.3 * y), 3 * y ./ (1 + 2 * y + 0.6 * x), 0, 2, 0, 3, 'Accuracy',30,'CutoffValues',[0 0], 'NullclineRange', [0 4 0 5])

% [alternatively:]

> f = @(x,y) 2 * x ./ (1 + x + 0.3 * y);

> g = @(x,y) 3 * y ./ (1 + 2 * y + 0.6 * x);

> augmented_implicit(f,g,0,2,0,3,'Accuracy',30,'CutoffValues',[0 0], 'NullclineRange', [0 4 0 5])
