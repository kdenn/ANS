function [rECI,vECI] = OEtoRVv2(e,i,Om,w,M,n,mu)
% Get ECI position and velocity from OE
% INPUT:
%   e: eccentricity
%   i: inclination [rad]
%   Om: right ascention [rad]
%   w: argument of perigee [rad]
%   M: mean anomaly [rad]
%   n: mean motion  [rad/s]

%% Get Ec
er = 100;
Eco = M;
iter = 0;
while er >= 1E-8 && iter <= 100
    del = -(Eco-e*sin(Eco)-M)/(1-e*cos(Eco));
    Ec = Eco + del;
    er = abs(Ec-Eco);
    Eco = Ec;
    iter = iter + 1;
end

%% Get PQW
a = (mu/n^2)^(1/3);
rPQW = [a*(cos(Ec)-e),a*sqrt(1-e^2)*sin(Ec),0]';
vPQW = (a*n)/(1-e*cos(Ec)).*[-sin(Ec),sqrt(1-e^2)*cos(Ec),0]';

%% Get Rotation
RzOm = [cos(-Om),sin(-Om),0;...
        -sin(-Om),cos(-Om),0;...
        0,0,1];
Rxi = [1,0,0;...
       0,cos(-i),sin(-i);...
       0,-sin(-i),cos(-i)];
Rzw = [cos(-w),sin(-w),0;...
        -sin(-w),cos(-w),0;...
        0,0,1];
R = RzOm*Rxi*Rzw;

%% Final ECI values
rECI = R*rPQW; % km
vECI = R*vPQW;

end