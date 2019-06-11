% TestSFM
% Kaitlin Dennison
% Stanford University
% Spring 2019
%{
    Test the structure from motion algorithm to obtain the 3D point of a
    crater that is matched in two images
%}

%% Setup
clear; clc;
TestConstants
sat1 = getGT(1,'TestScripts','eros3mill');
sat2 = getGT(2,'TestScripts','eros3mill');

%% Epipolar
% Camera Matrices
[R,t] = epiRot(sat1.rSatACAF,sat2.rSatACAF,sat1.rotACAFtoCam,sat2.rotACAFtoCam);
P = cell(1,2);
% P{1} = sat1.camFull;
% P{2} = sat2.camFull;
P{1} = sat1.camInt*[eye(3) zeros(3,1)];
P{2} = sat2.camInt*[R t];

% Image Sizes ([w;h])
imsize = flip([size(sat1.Im)', size(sat2.Im)']);

% Feature Locations ([c;r])
u = [1396 1328;
     1536 1587];
uc = u;
uc(1,:) = imsize(1,:) - uc(1,:); % flip the x-axis

%% Get correct ACI position
% Read in 3D mesh
if isequal(ast.filenameLR,'eros001708') || isequal(ast.filenameLR,'eros200700')
    [vACAF,f] = importTab([astModelFldr,ast.filenameLR,'.tab']);
    f = f+1;
else
    [vACAF,f] = importTab2([astModelFldr,ast.filenameLR,'.tab']);
end
vACI = rotVert(vACAF,ast.alp,ast.del,ast.W_0,ast.W_d,sat1.JD);
ftACI1 = transPt23(u(:,1)',sat1.rSatACI,sat1.rotACItoCam,f,vACI,sat1.camParams)'
points2D = transPt32(ftACI1',sat1.rSatACI,sat1.rotACItoCam,sat1.camParams)
ftACAF1 = sat1.rotACAFtoACI'*ftACI1;

%% SFM
% Xh = vgg_X_from_xP_nonlin(uc,P,imsize,[ftACAF1;1]); % ACAF frame
Xh = vgg_X_from_xP_lin(uc,P,imsize); % cam1 frame
rFeat = sat1.rotACAFtoACI*(sat1.rotACAFtoCam'*(Xh(1:3)./Xh(4))+sat1.rSatACAF)
uFeat = transPt32(rFeat',sat1.rSatACI,sat1.rotACItoCam,sat1.camParams)

%% Visual Check
figure()
imshow(sat1.Im)
hold on
plot(u(1,1),u(2,1),'rx')
plot(points2D(1),points2D(2),'c+')
plot(uFeat(1),uFeat(2),'go')
legend('Input','Ray Trace','SFM')
hold off







