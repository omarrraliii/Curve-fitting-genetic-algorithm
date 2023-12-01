# About the problem
Curve fitting is the process of constructing a
curve, or mathematical function (polynomial equation) 
that has the best fit to a series of data points, possibly subject to constraints.
In smooth curve fitting, the function is constructed to
approximately fit the data.
Given a set of points,
we would like to fit a curve to them using a polynomial equation.
# What the application do
This program is a genetic algorithm to find the best coefficients that would make the distance
between the polynomial function and the points minimum.
It takes the input from curve_fitting_input.txt in the following format:          
• First line: Number of datasets (must be at least 1)                                    
  For each dataset:                           
• Number of data points and polynomial degree separated by space                               
  For each point:                               
• x-coordinate and y-coordinate separated by space                                      
