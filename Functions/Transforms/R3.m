%{
Developed by Nathan Stacey
Last Updated - 12/20/2015

R3.m is a function that calculates the direction cosine matrix for
rotating a frame about the z-axis by some angle theta.
%}

function [R3] = R3(theta)
R3 = [cos(theta) sin(theta) 0;
      -sin(theta) cos(theta) 0;
      0 0 1]; 
end

