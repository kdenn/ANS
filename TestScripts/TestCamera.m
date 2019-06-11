% TestCamera
% Kaitlin Dennison
% Stanford University 
% Winter 2019
%{ 
    Test the camera projection to make sure the interinsics are correctly
    applied and the point of view of the satellite is correct.
%}

%% Setup
% Load simulation constants
TestConstants
sat = scANS;
cam = camOSIRIS;
[vACAF,f] = importTab([astModelFldr,ast.filename,'.tab']);

% Use the first index of time span
i = 1;

%% Dynamics
% Satellite position and velocity
[rSat,vSat] = sat.getRV(JDspan(i)); % ACI
R = rotWldCam(rSat,vSat,[0,0,0]); % Camera rotation matrix ACI -> cam

% Asteroid position and attitude
[rAst,vAst] = erosHCI(JDspan(i)); % HCI ecliptic
rAst = Reqec'*rAst; % SCI
vAst = Reqec'*vAst; % SCI
vACI = rotVert(vACAF,ast.alp,ast.del,ast.W_0,ast.W_d,JDspan(i));

%% Image
disp('Generating 2D image...')
tic;
Im = get2Dimage(f,vACI,rAst,rSat,R,cam);
t = toc;
disp([num2str(t),' seconds to generate'])

%% Validation
% Image size and orientation
if all(fliplr(size(Im)) == cam.pxNum)
    disp('Image size and orientation matches specs')
elseif all(size(Im) == cam.pxNum)
    disp('Image size matches specs but orientation is flipped')
else
    disp('Image size and orientation do not match specs')
end

% Perspective check
figure; hold on;
subplot(1,2,1); hold on
    imshow(Im)
hold off
subplot(1,2,2); hold on
    set(gca,'Color','k')
    patch('faces',f,'vertices',vACI,'FaceColor',0.6.*[1 1 1],'EdgeColor','none')
    material dull 
    light('Position',-rAst)
    axis([-20 20 -20 20 -20 20])
    view(rSat')
hold off














