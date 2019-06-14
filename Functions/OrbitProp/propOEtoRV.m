function [r,v] = propOEtoRV(oe,mu,t)
% Propagate from orbital elements to position and velocity in an inertial
% frame
% INPUT:
    % oe - orbital elements
        % [a e i Om w M_0]' 
        % (km - rad rad rad rad)
    % mu - gravitational parameter of the CB
    % t - time in seconds since epoch
% OUTPUT:
    % r - position (km) [3x1]
    % v - velocity (km/s) [3x1]
%-------------------------------------------------------------------------

%% Initialize
a = oe(1);
e = oe(2);
i = oe(3);
Om = oe(4);
w = oe(5);
n = sqrt(mu/(a^3));
M = mod(oe(6) + n*t, 2*pi);

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
r = R*rPQW; % km
v = R*vPQW;

end