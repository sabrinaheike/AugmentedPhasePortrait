
Version 1.0 

# AugmentedPhasePortrait
Matlab Code to create augmented phase portrait, introduced in S.H. Streipert, G.S.K. Wolkowicz: An augmented phase plane approach for discrete planar maps: Introducing next-iterate operators, https://doi.org/10.1016/j.mbs.2022.108924 [preprint: https://nam12.safelinks.protection.outlook.com/?url=http%3A%2F%2Farxiv.org%2Fabs%2F2210.07943&amp;data=05%7C01%7Csas887%40pitt.edu%7Cb1fdc08fa1a7492c59af08daafd5bb23%7C9ef9f489e0a04eeb87cc3a526112fd0d%7C1%7C0%7C638015629806079141%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C3000%7C%7C%7C&amp;sdata=HC6aRjnN8uMhCLHpqUollwkSvDdoda5nX9qL3wv0u0k%3D&amp;reserved=0], for planar difference equations x_(t+1)=f(x_t,y_t); y_(t+1)=g(x_t,y_t).


Input: (f,g,minx, maxx, miny, maxy, acc^*, cutoffx^*, cutoffy^*)
^* indicates optional parameters

f = functionhandle; right-hand side of the X-equation f = f(x,y)
g = functionhandle; right-hand side of the Y-equation g = g(x,y) 
minx = number; minimum x-value to be considered for plotting
maxx = number; maximum x-value to be considered for plotting
miny = number; minimum y-value to be considered for plotting
maxy = number; maximum y-value to be considered for plotting

optional parameters:
acc = number; optional parameter for number of arrows and signs of next-iterate operators; default value = 20
cutoffx  = number; if number given then a red dash-dotted curve is added in plot to show points (x,y) such that f(x,y)=cutoffx
cutoffy =  number; if number given then a orange dash-dotted curve is added in plot to show points (x,y) such that g(x,y)=cutoffy

Output: Plot of augmented phase portrait (dashed lines correspond to nullclines, solid lines to root-curves, signs to the sign of the next-iterate operator)


There are four programs, dependent on expressing the nullclines as functions in x or functions in y:
augmented_xx: expresses the nullclines of the X- and the Y-equation as functions in x.
augmented_xy: expresses the X-nullclines as functions in x and the Y-nullclines as functions in y. 
augmented_yx: expresses the X-nullclines as functions in Y and the Y-nullclines as functions in x.
augmented_yy: expresses the nullclines of the X- and the Y-equation as functions in y.

If a nullcline is expresses as a function in x, that is, y=l(x), then the corresponding next-iterate operator is L=g(x,y)-l(f(,x,y)).
If a nullcline is expresses as a function in y, that is, x=l(y), then the corresponding next-iterate operator is L=f(x,y)-l(g(,x,y)).

Note: Clear any x and y values you may have in your matlab memory by "clear x y" then type "syms x y" before calling the function

%%%%%%%%%%%%%%%%%%%%EXAMPLE 1: 

To get the augmented phase portrait for the system: 
x_(t+1)=2*x_t/(1+x_t+0.3*y_t),
y_(t+1)=3*y_t/(1+2*y_t+0.6*x_t)

plotted in [0,2]x[0,3] with default values, type:

> clear all

> syms x y

> augmented_xx(2*x/(1+x+0.3*y), 3*y/(1+2*y+0.6*x),0, 2,0, 3)

% [alternatively:]

> f=2*x/(1+x+0.3*y);

> g=3*y/(1+2*y+0.6*x);

> augmented_xx(f,g,0,2,0,3)

[For this example, the nullclines can be expressed as functions in x or y, so any of the other codes will also work: augmented_xy, augmented_yx, augmented_yy]

%%%%%%%%%%%%%%%%%%%%EXAMPLE 2: Fig. 11 in  S.H. Streipert, G.S.K. Wolkowicz: An augmented phase plane approach for discrete planar maps: Introducing next-iterate ...

% To get the augmented phase portrait for the system: 
 x_(t+1)=(1-0.4)*x_t-0.5*x_t*y_t+0.5,
 y_(t+1)=0.4*x_t*y_t
% plotted in [0,5]x[0,3] with fewer signs of the next-iterate operator and a line to check if solutions remain non-negative, type:

> clear all

> syms x y

> augmented_xy((1-0.4)*x-0.5*x*y+0.5, 0.5*x*y,0, 5,0, 3,15, 0,0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% For more examples, see Section 4 in S.H. Streipert, G.S.K. Wolkowicz: An augmented phase plane approach for discrete planar maps: Introducing next-iterate operators, submitted.

