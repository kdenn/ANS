% TestFeatDet
% Kaitlin Dennison
% Stanford University
% Winter 2019
%{
    Test the feature detection algorithm for obtaining feature positions
    and crater ellipses in 2D
%}

%% Setup
clear; clc;
TestConstants
sat = scANS;
cam = camOSIRIS;
GT = getGT(1,'TestScripts','eros3mill');

%% 2D Crater Detection
disp('Detecting craters...')
tic;
[craters,CC,ePar,Eog,Egr] = detectCraters(GT.Im,GT.rAstSCI,GT.rSatACI,GT.rotACItoCam,cam);
t = toc;
disp([num2str(t),' seconds to detect craters'])

%% Validation
disp([num2str(size(craters,1)),' craters detected'])
disp([num2str(CC.NumObjects),' edges detected'])

% Image with craters
figure; hold on
imshow(GT.Im); hold on
plotCraters(craters,CC,ePar,'r')
hold off
hold off

% Edges and grouping
figure; hold on
subplot(1,3,1); hold on
imshow(Eog)
title('Canny Edge Detection')
hold off
subplot(1,3,2); hold on
imshow(Egr)
title('Filtered Edge Detecion')
hold off
subplot(1,3,3); hold on
imshow(GT.Im)
for c = 1:CC.NumObjects
    if CC.LitSides(c)
        plot(CC.PixelSubList{c}(:,2),CC.PixelSubList{c}(:,1),'y')
    else
        plot(CC.PixelSubList{c}(:,2),CC.PixelSubList{c}(:,1),'c')
    end
end
title('Grouped Edges')
hold off
hold off








