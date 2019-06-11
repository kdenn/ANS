function [rHCI, vHCI] = erosHCI(JD)
% Returns Heliocentric position and velocity of eros on a given day
% INPUT:
%   JD: Julian date
% Output:
%   r: HCI position [x y z] [km]
%   v: HCI velocity [x' y' z'] [km/s]
% Notes:
%   - Checked by inputing JD of perihelion transit and checking rPQW within
%     the subfunction, OEtoRVv2.m
% -------------------------------------------------------------------------
% Kaitlin Dennison
% SLAB Research 
% Spring 2018
% -------------------------------------------------------------------------

%% Constants
mu_sun = 1.3271244004193938E11; % [km^3/s^2]
AU = 149597870.7; % AU conversion [km/AU]

% Orbital Elements (JPL HORIZONS)
Eros.JD_epoch = 2451170.5;                          % JD epoch [days]
Eros.e = 0.2228858603247133;                        % eccentricity
Eros.i = deg2rad(10.83015266864554);                % inclination [rad]
Eros.a = 1.458260038106518*AU;                      % semi-major axis [km]
Eros.w = deg2rad(178.6132327246327);                % arg or per [rad]
Eros.Om = deg2rad(304.4308844737856);               % RAAN [rad]
Eros.M_0 = deg2rad(208.1235381788443);              % mean anomaly [rad]
Eros.n = deg2rad(0.5596952381222570)/(24*60*60);    % mean motion [rad/s]
Eros.T = 2*pi/Eros.n;                               % period [s/rev]

% Rotation Parameters (IAU Report)
Eros.alp = 11.35;                                   % ascention [deg]
Eros.del = 17.22;                                   % declination [deg]
Eros.W_0 = 326.07;                                  % prime meridian [deg]
Eros.W_d = 1639.38864745;                           % rotation rate [deg]

%% Get r and v
t = (JD - Eros.JD_epoch) * (24*60*60); % given time in seconds since epoch
M = mod(Eros.M_0 + Eros.n * t, 2 * pi);
[rHCI,vHCI] = OEtoRVv2(Eros.e,Eros.i,Eros.Om,Eros.w,M,Eros.n,mu_sun);
% rHCI = rotEQUtoECL()*rHCI;
% vHCI = rotEQUtoECL()*vHCI;

end
