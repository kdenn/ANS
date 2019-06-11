function v_rot = rotVert(v,alp,del,W_0,W_d,JD,rev)
% Rotate the vertices from ACAF -> ACI frame
% INPUT:
%   v: matrix of vertices [nx3] [x y z] in Equatorial Plane
%   alp: alpha, right ascension [deg]
%   del: delta, declination [deg]
%   W_0: prime meridian [deg]
%   W_d: rotation rate [deg/day]
%   JD: Julian Date
% OUTPUT:
%   v_rot: matrix of vertices [nx3] [x y z] in Ecliptic Plane
% NOTES:
%   Standard epoch: J2000.0 = JD 2451545.0 (2000 Jan 1 12:00)
% -------------------------------------------------------------------------
% Kaitlin Dennison [5/1/2018]
% Stanford University - Space Rendezvous Laboratory
% -------------------------------------------------------------------------
if nargin == 6
    rev = false;
end
d = JD - 2451545.0; % days since standard epoch
n = size(v,1);
v_rot = zeros(n,3);
Z1 = R3(deg2rad(-W_0-W_d*d));
X1 = R1(deg2rad(del-90));
Z2 = R3(deg2rad(-alp-90));
R = Z2*X1*Z1;
if rev
    for i = 1:n
        v_rot(i,:) = (R'*rotEQUtoECL()'*v(i,:)')';
    end
else
    for i = 1:n
        v_rot(i,:) = (rotEQUtoECL()*R*v(i,:)')';
    end
end
end