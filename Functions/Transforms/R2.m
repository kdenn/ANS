%{
Developed by Nathan Stacey
Last Updated - 12/20/2015

R2.m is a function that calculates the direction cosine matrix for
rotating a frame about the y-axis by some angle theta.
%}

function [R2] = R2(theta)
R2 = [cos(theta) 0 -sin(theta);
      0 1 0;
      sin(theta) 0 cos(theta)];
end

