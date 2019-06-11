function R = rotWldCam(rSat,vSat,aSat)
% Rotation matrix from ACI to Cam
h = cross(rSat,vSat); % angular momentum of sat
% AA273 AR lec 5 s10
kc = rSat./norm(rSat); % Directly towards asteroid CM
jc = h./norm(h); % Aligned with angular momentum
ic = cross(jc,kc)./norm(cross(jc,kc)); % Complete the RH triad
% Rotate from looking directly at CM via 3-2-1
Rc = R1(aSat(1))*R2(aSat(2))*R3(aSat(3));
iu = Rc*ic;
ju = Rc*jc;
ku = Rc*kc;
iw = [1 0 0];
jw = [0 1 0];
kw = [0 0 1];
R = [dot(iw,iu) dot(jw,iu) dot(kw,iu);
     dot(iw,ju) dot(jw,ju) dot(kw,ju);
     dot(iw,ku) dot(jw,ku) dot(kw,ku)];
end