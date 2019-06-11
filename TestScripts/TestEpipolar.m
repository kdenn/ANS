% TestEpipolar
% Kaitlin Dennison
% Stanford University
% Spring 2019
%{
    Test the epipolar matching algorithm to match feature in sequential
    images
%}

%% Setup
clear; clc;
TestConstants
sat1 = getGT(1,'TestScripts','eros3mill');
sat2 = getGT(2,'TestScripts','eros3mill');

%% Detect Features
% [craters1,CC1,ePar1,Eog1,Egr1] = detectCraters(Im1,rAst1,rSat1,rotWldCam(rSat1,vSat1,aSat1),camParams1);
% [craters2,CC2,ePar2,Eog2,Egr2] = detectCraters(Im2,rAst2,rSat2,rotWldCam(rSat2,vSat2,aSat2),camParams2);
load('TestScripts\craters.mat')
u = [1396 1328;
     1536 1587];

%% Epipolar Lines
imsize = flip([size(sat1.Im)', size(sat2.Im)']);
P = cell(1,2);
P{1} = sat1.camFull;
P{2} = sat2.camFull;
F = vgg_F_from_P(P);
draw_epipolar(u(1,1),u(2,1),F,sat1.Im,sat2.Im)

%% Visual Check
% figure; hold on
% subplot(1,2,1); hold on
%     imshow(Im1); hold on
%     for c = 1:size(craters1,1)
%         ellipse(craters1(c,3),craters1(c,4),craters1(c,5),craters1(c,1),craters1(c,2))
%     end
%     hold off
%     hold off
% subplot(1,2,2); hold on
%     imshow(Im2); hold on
%     for c = 1:size(craters2,1)
%         ellipse(craters2(c,3),craters2(c,4),craters2(c,5),craters2(c,1),craters2(c,2))
%     end
%     hold off
%     hold off





