function [rECI,vECI] = OEtoRV(a,e,i,Om,w,nu,mu)
% Get the ECI r and v from OEs
% INPUT:
%   a: semi-major axis [km]
%   e: eccentricity
%   i: inclination [rad]
%   Om: RAAN [rad]
%   w: arg of periapsis [rad]
%   nu: true anomaly [rad]
%   mu: gravitational parameter of central body
% -------------------------------------------------------------------------

Ec = acos((e+cos(nu))/(1+e*cos(nu))); 
if nu > pi
    Ec = 2*pi-Ec;
end
n = sqrt(mu/a^3);
% Perifocal Coordinates
rPQW = [a*(cos(Ec)-e);
        a*sqrt(1-e^2)*sin(Ec);
        0];
vPQW = (a*n)/(1-e*cos(Ec)).*[-sin(Ec);
                             sqrt(1-e^2)*cos(Ec);
                             0];
% Rotation Matrices
RzOm = [cos(-Om),sin(-Om),0;...
        -sin(-Om),cos(-Om),0;...
        0,0,1];
Rxi = [1,0,0;...
       0,cos(-i),sin(-i);...
       0,-sin(-i),cos(-i)];
Rzw = [cos(-w),sin(-w),0;...
        -sin(-w),cos(-w),0;...
        0,0,1];
if i == 0 && e ~= 0
    % Equatorial & Elliptical
    R = Rzw;
elseif e == 0 && i ~= 0
    % Circular & Inclined
    R = RzOm*Rxi;
elseif i == 0 && e == 0
    % Equatorial & Circular
    R = 1;
else
    R = RzOm*Rxi*Rzw;
end
% Rotate to ECI
rECI = R*rPQW; % km
vECI = R*vPQW;
end


