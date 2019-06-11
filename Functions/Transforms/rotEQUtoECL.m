function R = rotEQUtoECL()
% Rotation matrix to convert from equatorial to ecliptice plane
tilt = 23.43688; % [deg]
R = [1 0           0;
     0 cosd(tilt)  sind(tilt);
     0 -sind(tilt) cosd(tilt)];