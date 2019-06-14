function eul = rotMtoEul(R)
% Get the Euler angles from the 3x3 rotation matrix using 3-1-3 convention
% R = R3(a1)*R1(a2)*R3(a3)

a2 = acos(R(3,3));
a1 = acos(R(2,3)/sin(a2));
a3 = asin(R(3,1)/sin(a2));

eul = [a1;a2;a3];

end