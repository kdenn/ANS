function GT = getGT(num,fldr,model)

Im = imread([fldr,'/',model,'_im',num2str(num),'.png']);
meta = csvread([fldr,'/',model,'_da',num2str(num),'.csv']);
JD = datenum(meta(4,3),meta(4,4),meta(4,5),0,0,meta(4,6))+1721058.5;

% Camera Parameters
camParams = meta(3,:);
imPxAct = meta(5,[4,5]);
camParams(7:8) = imPxAct;
camInt = getCamIntrinsic(meta(3,:));

% Satellite and Asteroid positions
rSat = meta(1,1:3)'; % ACI
vSat = meta(1,4:6)'; % ACI
rAst = meta(2,1:3)'; % SCI
vAst = meta(2,4:6)'; % SCI

% ACI ecliptic to camera frame rotation matrix
aSat = meta(5,1:3)'; % rad off from ast center
R = rotWldCam(rSat,vSat,aSat);
RFI = rotACAFtoACI(JD);
P = camInt*[R*RFI R*RFI*(-(RFI'*rSat))];

% Put into a structure
GT.Im = Im; % image
GT.meta = meta; % meta data matrix
GT.JD = JD; % Juldian date of image
GT.camParams = camParams; % camera parameters vector
GT.camInt = camInt; % camera intrinsic matrix
GT.camFull = P; % camera full matrix to relate ACAF to Cam
GT.rSatACI = rSat; % sat ACI position
GT.rSatACAF = RFI'*rSat;  % sat ACAF position
GT.rAstSCI = rAst; % asteroid SCI position
GT.vAstSCI = vAst; % asteroid SCI velocity
GT.rotACItoCam = R; % rot matrix ACI to camera
GT.rotACAFtoACI = RFI; % rot matrix ACAF to ACI
GT.rotACAFtoCam = R*RFI; % rot matrix from ACAF to Cam

end