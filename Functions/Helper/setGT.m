function GT = setGT(sat,ast,JD,Im)



%% 
if nargin == 3
    Im = [];
end
    
%% Satellite Position
M = mod(sat.M_0 + sat.n*(JD-sat.JD_epoch)*60*60*24,2*pi);
Ec = iterMtoEc(M,sat.e,1E-8);
nu = acos((cos(Ec)-sat.e)/(1-sat.e* cos(Ec)));
if Ec > pi
    nu = 2*pi-nu;
end
[rSat,vSat] = OEtoRV(sat.a,sat.e,sat.i,sat.Om,sat.w,nu,ast.mu); % ACI

%% Generate Images
[rAst,vAst] = erosHCI(JDspan(i)); % HCI ecliptic
rAst = Reqec'*rAst; % SCI
vAst = Reqec'*vAst; % SCI
vACI = rotVert(vACAF,ast.alp,ast.del,ast.W_0,ast.W_d,JDspan(i));
R = rotWldCam(rSat,vSat,[0,0,0]);
% Create a csv file with the ground truth data that can be read
% in in conjunction with the images.
%         if dispStatus, disp('    Saving Meta-Data'); end
[Ye,Mo,Da,Hr,Mi,Se] = datevec(JDspan(i)-1721058.5);
Se = Se + (Hr*60*60) + (Mi*60);
% Lighting Direction
dat = transPt32(-rAst',rSat,R,camParams)';
Ldxn = -dat'./norm(dat);
a = acosd(dot(rSat,-rAst)/(norm(rSat)*norm(-rAst)));
if a > 90
    Ldxn = -Ldxn;
end
% Write file
meta = [rSat', vSat';
    rAst', vAst';
    camParams(1:6);
    Ldxn, Ye, Mo, Da, Se;
    0, 0, 0, size(Im), 0];


%% Put into a structure
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