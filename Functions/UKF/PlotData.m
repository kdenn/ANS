% PlotData
% Kaitlin Dennison
% Stanford University
% Winter 2019
%{
    Build a database of images and features
%}

%% Setup
clear; clc;

% Load simulation constants
SimConstants

% Load simulation data
load([saveToFldr,saveWorkspace]); % Set of Features
load([saveToFldr,ftDBloc]); % Set of Features [x,y,z,JD,count,Ptri,Mx,My,sum_a,sum_b,sum_alpha]
rmv = ftDB(:,5) < 3;
% ftDB(rmv,:) = [];
nf = size(ftDB,1);
ni = 400;

% Load Mesh
if isequal(ast.filenameLR,'eros001708') || isequal(ast.filenameLR,'eros200700')
    [vACAF,f] = importTab([astModelFldr,ast.filenameLR,'.tab']);
    f = f+1;
else
    [vACAF,f] = importTab2([astModelFldr,ast.filenameLR,'.tab']);
end

%% 
x_sum = zeros(ni,3); % x_ACAF
P_sum = zeros(ni,3); % x_ACAF 3D Covariance
m_sum = zeros(ni,1); % Mahalanobis Distance
r_sum = zeros(ni,1); % Residuals
M_sum = zeros(ni,2); % z_px 2D Covariance
c_sum = zeros(100,2); % z_px 2D Covariance per correlation count
ct = zeros(ni,3); % correlation count per image
ct_c = zeros(100,1); % correlation count

for ftNum = 1:nf
    dat = csvread([saveToFldr,'FeatureInfo/',ast.filename,'_ft',num2str(ftNum),'.csv']);
    % dat = [ft#,c,r,a,b,alpha,score,mah,xACI,xACAF,diag(P)',diag(M)'];
    x_f = dat(end,13:15);
%     c_sum(1,:) = c_sum(1,:) + dat(1,19:20);
%     ct_c(1) = ct_c(1) + 1;
    if size(dat,1) < 2 || size(dat,2) < 20
        continue
    end
    for c = 2:size(dat,1)
        im = dat(c,1);
        x_sum(im,:) = x_sum(im,:) + x_f - dat(c,13:15);
        P_sum(im,:) = P_sum(im,:) + sqrt(dat(c,16:18));
        M_sum(im,:) = M_sum(im,:) + sqrt(dat(c,19:20));
        m_sum(im,:) = m_sum(im,:) + dat(c,8);
        r_sum(im,:) = r_sum(im,:) + dat(c,9);
        c_sum(c,:) = c_sum(c,:) + sqrt(dat(c,19:20));
        ct(im,:) = ct(im,:) + [1,1,1];
        ct_c(c) = ct_c(c) + 1;
    end
end

sz = camParams(2)*camParams(3);
i_avg = zeros(ni,1);
in_ang = zeros(ni,1);
for i = 1:ni
    Im = im2double(imread([saveToFldr,ast.filename,'_im',num2str(i),'.png']));
    i_avg(i) = sum(sum(Im))/sz;
    dat = csvread([saveToFldr,ast.filename,'_da',num2str(i),'.csv']);
    rSat = dat(1,1:3)';
    rSun = -dat(2,1:3)';
    in_ang(i) = mod(acosd(dot(rSat,rSun)/(norm(rSat)*norm(rSun))),360);
end

%% Robust Crater
[m,ftBest] = max(ftDB(:,5));
dat = csvread([saveToFldr,'FeatureInfo/',ast.filename,'_ft',num2str(ftNum),'.csv']);
% dat = [im#,c,r,a,b,alpha,score,mah,xACI,xACAF,diag(P)',diag(M)'];
x_f = dat(end,13:15);

figure; hold on
    subplot(3,1,1); hold on
        plot(dat(:,1),1000.*(dat(:,13)-x_f(1)),'r+')
        plot(dat(:,1),3000.*sqrt(dat(:,16)),'b')
        plot(dat(:,1),-3000.*sqrt(dat(:,16)),'b')
        ylabel('X_{ACAF} residual [m]')
        axis([1,dat(end,1),-500,500])
        hold off
    subplot(3,1,2); hold on
        plot(dat(:,1),1000.*(dat(:,14)-x_f(2)),'r+')
        plot(dat(:,1),3000.*sqrt(dat(:,17)),'b')
        plot(dat(:,1),-3000.*sqrt(dat(:,17)),'b')
        ylabel('Y_{ACAF} residual [m]')
        axis([1,dat(end,1),-500,500])
        hold off
    subplot(3,1,3); hold on
        plot(dat(:,1),1000.*(dat(:,15)-x_f(3)),'r+')
        plot(dat(:,1),3000.*sqrt(dat(:,18)),'b')
        plot(dat(:,1),-3000.*sqrt(dat(:,18)),'b')
        ylabel('Z_{ACAF} residual [m]')
        axis([1,dat(end,1),-500,500])
        xlabel('Image Number')
        hold off
    hold off
hold off

%% Covariance
x_avg = 1000.*x_sum./ct;
P_avg = 1000.*P_sum./ct;
M_avg = M_sum./ct(:,2);
m_avg = m_sum./ct(:,1);
r_avg = r_sum./ct(:,1);
idx = (1:ni)';
idx(ct(:,1)==0,:) = [];
x_avg(ct(:,1)==0,:) = [];
P_avg(ct(:,1)==0,:) = [];
m_avg(ct(:,1)==0,:) = [];
r_avg(ct(:,1)==0,:) = [];
M_avg(ct(:,1)==0,:) = [];

%{
figure; hold on
    subplot(3,1,1); hold on
        plot(idx,x_avg(:,1),'r')
        plot(idx,3.*P_avg(:,1),'b')
        plot(idx,-3.*P_avg(:,1),'b')
        ylabel('X_{ACAF} residual [m]')
        hold off
    subplot(3,1,2); hold on
        plot(idx,x_avg(:,2),'r')
        plot(idx,3.*P_avg(:,2),'b')
        plot(idx,-3.*P_avg(:,2),'b')
        ylabel('Y_{ACAF} residual [m]')
        hold off
    subplot(3,1,3); hold on
        plot(idx,x_avg(:,3),'r')
        plot(idx,3.*P_avg(:,3),'b')
        plot(idx,-3.*P_avg(:,3),'b')
        ylabel('Z_{ACAF} residual [m]')
        xlabel('Image Number')
        hold off
    hold off
hold off
%}

%% Mahalanobis

%{
figure; hold on
    subplot(3,1,1); hold on
        plot(idx,3.*M_avg(:,1),'b')
        plot(idx,-3.*M_avg(:,1),'b')
        ylabel('X_{2D} [px]')
        hold off
    subplot(3,1,2); hold on
        plot(idx,3.*M_avg(:,2),'b')
        plot(idx,-3.*M_avg(:,2),'b')
        ylabel('Y_{2D} [px]')
        hold off
    subplot(3,1,3); hold on
        plot(idx,m_avg,'r')
        ylabel('Mahalanobis')
        xlabel('Image Number')
        hold off
%     subplot(4,1,4); hold on
%         plot(idx,r_avg,'r')
%         ylabel('Residual')
%         hold off
hold off
%}

%% Correlation Counts
c_avg = c_sum./ct_c;
figure; hold on
subplot(3,1,1); hold on
    histogram(ftDB(:,5))
    title('Correlations per Feature')
    xlabel('Number of Correlations')
    ylabel('Number of Features')
    hold off
subplot(3,1,2); hold on
    plot((1:length(ct_c))',c_avg)
    title('Average 2D Cov per Seq Correlation')
    xlabel('Number of Correlations')
    ylabel('\Sigma_{2D} [px]')
    legend('\Sigma_x','\Sigma_y')
    hold off
subplot(3,1,3); hold on
    title('Correlations per Image')
    xlabel('Image Number')
    yyaxis right
    plot((1:ni)',180.*i_avg./max(i_avg))
    plot((1:ni)',in_ang)
    ylabel('Incidence Ang')
    yyaxis left
    plot((1:ni)',ct(1:ni))
    ylabel('Corr per Im')
    hold off
hold off

%%
figure; hold on
    subplot(3,1,1); hold on
        scatter(in_ang,ct(1:ni),'*')
        xlabel('Incidence Angle')
        ylabel('Correlations')
        hold off
    subplot(3,1,2); hold on
        scatter(i_avg,ct(1:ni),'*')
        xlabel('Average Image Intensity')
        ylabel('Correlations')
        hold off
    subplot(3,1,3); hold off
        scatter(in_ang,i_avg,'*','CData',ct(1:ni))
        xlabel('Incidence Angle')
        ylabel('Average Image Intensity')
        colorbar
        hold off
    hold off

%% 3D Features


figure; hold on
    set(gca,'Color','k')
    patch('faces',f,'vertices',vACAF,'FaceColor',0.6.*[1 1 1],'EdgeColor','none')
    material dull 
    light('Position',[1000000000,0,0])
    axis([-20 20 -20 20 -20 20])
    axis equal
    scatter3(ftDB(:,1),ftDB(:,2),ftDB(:,3),'*','CData',ftDB(:,5))
    colorbar
    view(3)
%     plot3(ftDB(ftMx,1),ftDB(ftMx,2),ftDB(ftMx,3),'g*')
hold off