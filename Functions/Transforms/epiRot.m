function [R,t] = epiRot(rSat1f,rSat2f,Rot1,Rot2)
% Rotation matrix and position vector from sat1 to sat2
% INPUT:
%   rSat#f - ACAF position of sat#
%   Rot# - 3x3 rotation matrix from ACAF to cam#
% OUTPUT:
%   R - 3x3 rotation matrix from cam1 to cam2
%   t - 3x1 position vector from cam2 to cam1 in the cam2 frame

tACAF = rSat1f - rSat2f;
t = Rot2*tACAF;
R = Rot2*Rot1';