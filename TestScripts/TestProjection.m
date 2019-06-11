% TestProjection
% Kaitlin Dennison
% Stanford University
% Spring 2019
%{
    Test the projections of feature centroids between 2D and 3D
%}

%% Setup
%{
clear; clc;
TestConstants
sat = scANS;
cam = camOSIRIS;
GT = getGT(1,'TestScripts','eros3mill');
load('TestScripts\craters.mat')
[vACAF,f] = importTab([astModelFldr,ast.filename,'.tab']);
vACI = rotVert(vACAF,ast.alp,ast.del,ast.W_0,ast.W_d,GT.JD);
disp('Loaded')
%}

%% 2D Crater Detection
disp('2D to 3D...')
tic
craters2 = transPt23(craters1,GT.rSatACI,GT.rotACItoCam,f,vACI,cam);
t = toc;
disp([num2str(t), ' seconds to project'])
craters3 = transPt32(craters2,GT.rSatACI,GT.rotACItoCam,cam);

%% Validation

% Edges and grouping
figure; hold on
    subplot(1,2,1); hold on
        imshow(GT.Im); hold on
            scatter(craters1(:,1),craters1(:,2),'*r')
            scatter(craters3(:,1),craters3(:,2),'*b')
            plotCraters(craters1)
            title('Original 2D')
        hold off
    hold off
    subplot(1,2,2); hold on
        set(gca,'Color','k')
        patch('faces',f,'vertices',vACI,'FaceColor',0.6.*[1 1 1],'EdgeColor','none')
        material dull 
        light('Position',-GT.rAstSCI)
        axis([-20 20 -20 20 -20 20])
        view(GT.rSatACI')
        scatter3(craters2(:,1),craters2(:,2),craters2(:,3),'*r')
    hold off
hold off








