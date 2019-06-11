%{
Developed by Nathan Stacey
Last Updated - 12/20/2015

R1.m is a function that calculates the direction cosine matrix for
rotating a frame about the x-axis by some angle theta.
%}

function [R1] = R1(theta)
R1 = [1 0 0;
      0 cos(theta) sin(theta);
      0 -sin(theta) cos(theta)];
end

