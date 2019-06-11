function R = rotACAFtoACI(JD)
% Rotatation matrix from ACAF -> ACI

%% Eros Constants
alp = 11.35;                                   % ascention [deg]
del = 17.22;                                   % declination [deg]
W_0 = 326.07;                                  % prime meridian [deg]
W_d = 1639.38864745;                           % rotation rate [deg]

%% Calculations
d = JD - 2451545.0; % days since standard epoch
Z1 = R3(deg2rad(-W_0-W_d*d));
X1 = R1(deg2rad(del-90));
Z2 = R3(deg2rad(-alp-90));
R = Z2*X1*Z1;

end