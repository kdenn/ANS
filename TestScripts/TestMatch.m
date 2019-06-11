% TestMatch
% Kaitlin Dennison
% Stanford University
% Spring 2019
%{
    Test the matching algorithms
%}

%% Setup
clear; clc;
TestConstants
sat = scANS;
cam = camOSIRIS;
load(ftDBloc)
load('TestScripts\craters.mat')
pt1 = getGT(1,'TestScripts','eros3mill');
pt2 = getGT(2,'TestScripts','eros3mill');

%% Mahalanobis
% Picture 1
ftACI1 = rotVert(ftDB(:,1:3),ast.alp,ast.del,ast.W_0,ast.W_d,pt1.JD);
ftDB1 = ftDB;
for k = 1:size(ftDB1,1)
    % M = get2Dcov(ftACI(k,1:3),R,ftDB(k,[6,9,11])',rSat,camParams);
    P0 = getDBfeat(ftDB1,k,'3Dcov');
    M = get2Dcov(ftACI1(k,1:3),P0,@(xs) hmeas(xs,pt1.JD,ast,pt1.rSatACI,pt1.rotACItoCam,cam));
    if M(1,1) < tolM
        M = (sqrt(M) + eye(2).*40).^2;
    end
    ftDB1 = setDBfeat(ftDB1,k,'2Dcov',M);
end
%     smlM = ftDB(:,12) < tolM;
%     ftDBinfl(smlM,12:13) = ones(sum(smlM),2).*(40^2); % inflate M if too small
[corr1,uncorr1] = corrFeat(craters1,pt1.rSatACI,pt1.rotACItoCam,ftACI1,ftDB1,cam,sigM,maxR); % [ft#,c,r]

% Picture 2
ftACI2 = rotVert(ftDB(:,1:3),ast.alp,ast.del,ast.W_0,ast.W_d,pt2.JD);
ftDB2 = ftDB;
for k = 1:size(ftDB2,1)
    % M = get2Dcov(ftACI(k,1:3),R,ftDB(k,[6,9,11])',rSat,camParams);
    P0 = getDBfeat(ftDB2,k,'3Dcov');
    M = get2Dcov(ftACI2(k,1:3),P0,@(xs) hmeas(xs,pt2.JD,ast,pt2.rSatACI,pt2.rotACItoCam,cam));
    if M(1,1) < tolM
        M = (sqrt(M) + eye(2).*40).^2;
    end
    ftDB2 = setDBfeat(ftDB2,k,'2Dcov',M);
end
%     smlM = ftDB(:,12) < tolM;
%     ftDBinfl(smlM,12:13) = ones(sum(smlM),2).*(40^2); % inflate M if too small
[corr2,uncorr2] = corrFeat(craters2,pt2.rSatACI,pt2.rotACItoCam,ftACI2,ftDB2,cam,sigM,maxR); % [ft#,c,r]

%% Visual Check
ftDB2D1 = transPt32(ftACI1,pt1.rSatACI,pt1.rotACItoCam,pt1.camParams);
plotMatch(pt1.Im,craters1,ftDB2D1,ftDB1,corr1,1)
ftDB2D2 = transPt32(ftACI2,pt2.rSatACI,pt2.rotACItoCam,pt2.camParams);
plotMatch(pt2.Im,craters2,ftDB2D2,ftDB2,corr2,1)