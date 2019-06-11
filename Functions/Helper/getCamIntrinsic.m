function K = getCamIntrinsic(camParams)
% Retreive the camera intrinsic matrix from the camParams vector

%% Camera Constants
foc = camParams(1); % [mm] focal length
pxNum = camParams(2:3); % pixels in CCD pixel array rxc
pxSze = camParams(4:5); % [mm] wxh of one px
imSze = pxSze.*pxNum([2,1]); % [mm] size of CCD array
k = 1./pxSze; % [px/mm]

%% Get Camera Intrinsics
% AA273 AR lec 5
u0 = imSze(1)/2;
v0 = imSze(2)/2;
K = [foc*k(1) 0        k(1)*u0;
     0        foc*k(2) k(2)*v0;
     0        0        1];
 
end