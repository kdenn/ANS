function points2D = transPt32(points3D,rSat,R,camera)
% Convert 3D ACI crater positions to 2D image pixel frame
% INPUT:
    % points3D: [x y z] positions in ACI ecliptic frame [nx3]
    % rSat: position of satellite in ACI Ecliptic [3x1]
    % R: 3x3 rotation matrix from ACI Ecliptic to Camera (see rotWldCam)
    % camParams: array of camera parameters
        % [foc, pxNum, pxSze, radDist, imPxAct]
% OUTPUT:
    % points2D: pixel col/row coordinates in image [c r] [nx2]

%% Parse Input
nft = size(points3D,1);
points2D = zeros(nft,2);

%% Get Camera Intrinsics
% AA273 AR lec 5
pxNum = camera.pxNum; % pixels in CCD pixel array rxc
t = R*(-rSat);
P = camera.A*[R t];

%% Convert from World to 2D Camera
% AA273 AR lec 5 s13
for i = 1:nft
    Pw = points3D(i,1:3)';
    Pwh = [Pw; 1];
    ph = P*Pwh;
    p = ph(1:2)./ph(3);
    points2D(i,1:2) = p';
end
% Put into row-col frame
points2D(:,1) = pxNum(1) - points2D(:,1);
end

















