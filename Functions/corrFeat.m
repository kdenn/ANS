function [corr,uncorr] = corrFeat(craters_new,rSat,R,ftACI,ftDB,camera,sigM,maxR)
% Correlate detected features to the current database
% INPUT:
% craters: matrix with detected crater data [nx6] [c r a b alpha score]
% rSat: position of satellite in ACI Equatorial [3x1]
% R: 3x3 rotation matrix from ACI Equatorial to Camera (see rotWldCam)
% camera: camera object (see camera.m class)
% [foc, pxNum, pxSze, radDist, imPxAct]
% f: faces defined by vertex numbers [nfx3] [v1 v2 v3]
% vACI: matrix of vertices [nvx3] [x y z] ACI Equatorial
% ftACI: feature database [nftx3] [x y z] ACI Equatorial
% dr: radius of search for matching features
% M: Mahalanobis distance matrix [3x3], essentially a covariance matrix
% in the ACI frame. You could probably replace this with the 3x3
% covariance matrix for the ftACI position from the UKF.
% OUTPUT:
% corr: correlated features and their px coordinates [ft# c r a b alpha score,mah]
% uncorr: uncorrelated features with new feature numbers [ft# c r a b alpha score]
% ftACI_new: ACI form of the new features (empty if using 2D)

%% Parse Input
nc = size(craters_new,1);
nf = size(ftACI,1);
if nf == 0
    corr = zeros(0,9);
    uncorr = [(1:nc)', craters_new];
    return
end
C = 2;
D = 15*pi/180;

%% What features are visible?
vizFaces = zeros(nf,1);
for i = 1:nf
    rFeatSat = rSat-ftACI(i,:)';
    incidence = acos(dot(ftACI(i,:)./norm(ftACI(i,:)),rFeatSat./norm(rFeatSat)));
    if incidence > pi/2
        vizFaces(i) = 0;
    else
        vizFaces(i) = 1;
    end
end
vizMap = (1:(length(vizFaces)))';
vizMap(~vizFaces) = [];
ftACI_viz = ftACI(vizMap,:);
ftDB_viz = ftDB(vizMap,:);
nfv = length(vizMap);

%% Setup
corr = zeros(nfv,9); % [ft# c r a b alpha score mah res]
corr(:,1) = (1:nfv)';
corr(:,[8,9]) = ones(nfv,2).*inf;
uncorr = zeros(0,7);
craters = transPt32(ftACI_viz,rSat,R,camera);
[X,Y] = meshgrid(craters(:,1),craters(:,2));
dist = sqrt((X-X').^2+(Y-Y').^2);
for n = 1:nc
    z = craters_new(n,1:2)';
    Mah_best = inf;
    m_best = 0;
    for m = 1:nfv
        d = z-craters(m,1:2)';
        cov2D = [ftDB_viz(m,12:13);ftDB_viz(m,13:14)]; % + eye(2).*craters_new(n,6)
        Mah = sqrt(d'*(cov2D\d));
        if Mah > sigM
            continue
        end
        if Mah < Mah_best
            Mah_best = Mah;
            m_best = m;
        end
    end
    if m_best == 0
        uncorr = [uncorr;
            nf+1,craters_new(n,:)]; % [ft#,c,r,a,b,alpha,score]
        nf = nf + 1;
    else
        % Check residuals to other features
        dn = sqrt((z(1)-craters(:,1)).^2+(z(2)-craters(:,2)).^2);
        res = norm(dist(:,m_best) - dn);
        if (Mah_best < corr(m_best,8)) && res < corr(m_best,9) && res < maxR
            corr(m_best,:) = [m_best,craters_new(n,:),Mah_best,res]; % [ft#,c,r,a,b,alpha,score]
        end
    end
end

corr(:,1) = vizMap;
corr((corr(:,9)>=maxR),:) = [];

end