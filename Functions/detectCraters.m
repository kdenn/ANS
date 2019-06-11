function [craters,CC,ePar,Eog,Egr,crit] = detectCraters(Im,rAst,rSat,R,camera,Params)
% Detect the craters in an image of an asteroid.
% INPUT:
    % Im: grayscale image
    % rAst: position of satellite in SCI Equatorial [3x1]
    % rSat: position of satellite in ACI Equatorial [3x1]
    % R: 3x3 rotation matrix from ACI Equatorial to Camera (see rotWldCam)
    % camParams: array of camera parameters
        % [foc, pxNum, pxSze, radDist, imPxAct]
    % Params: (optional) structure with detection parameters, if you leave
    % this out of supply a structure with only a few of the parameters, the
    % defaults will be filled in. Defaults:
        % dThreshUp = 3; % distance between centers
        % dThreshLo = 0.01;
        % ceSml = 0; % sobel edge detection
        % ceLrg = 0.02;
        % lthresh = 4; % length ratio
        % sThresh = 0.02; % slope of linear fit
        % gThresh = 0.64; % angle from lighting direction
        % rThresh = 10; % norm of residuals of linear fit
        % wr = 1; % weight on length ratio wr*rto
        % ws = 1; % weight on slope ws/s
        % wd = 5; % weight on distance wd*d
        % wn = 0; % weight on norm of residuals wn*n
% OUTPUT:
    % craters: matrix with detected crater data [nx6] [c r a b alpha score]
    % CC: connected component structure with all of the edge data
    % associated with the craters
        % Connectivity: 8
        % ImageSize: [numRow numCol]
        % NumObjects: Ne, number of edges detected
        % PixelIdxList: {1×Ne cell}
        % PixelSubList: {1×Ne cell}
        % DistToLight: {1×Ne cell}
        % CenterPx: [1×Ne double]
        % ClosestPx: [1×Ne double]
        % FarthestPx: [1×Ne double]
        % LitSides: [1×Ne double]
        % MatchedEdges: [Ne×1 double]
    % ePar: edge pair matching info [nC×7 double]
        % [e1 e2 rto s normr d cst];
    % Eog: binary mask of edges before gradient mask applied
    % Gck: binary mask with gradient direction thresholding
    % Egr: binary mask of edges with Gck applied

if nargin == 5
    % Default Parameters
    Params.lthresh = 2.5; % length ratio
    Params.dThreshUp = 3.7; % distance between centers
    Params.dThreshLo = 0.01;
    Params.ceSml = 0.2; % sobel edge detection
    Params.ceLrg = 0.4;
    Params.sThresh = 6.5; % slope of linear fit
    Params.gThresh = 0.64; % angle from lighting direction
    Params.rThresh = 10; % norm of residuals of linear fit
    Params.wr = 1; % weight on length ratio wr*rto
    Params.ws = 1; % weight on slope ws/s
    Params.wd = 5; % weight on distance wd*d
    Params.wn = 0; % weight on norm of residuals wn*n
    Params.smoothing = true; % turn on smoothing
    Params.gaussSig = 7; %5/(norm(rSat)/35); % sigma for Gaussian blur
else
    if ~isfield(Params,'lthresh'), Params.lthresh = 4; end
    if ~isfield(Params,'dThreshUp'), Params.dThreshUp = 3; end
    if ~isfield(Params,'dThreshLo'), Params.dThreshLo = 0.01; end
    if ~isfield(Params,'ceSml'), Params.ceSml = 0.2; end
    if ~isfield(Params,'ceLrg'), Params.ceLrg = 0.4; end
    if ~isfield(Params,'sThresh'), Params.sThresh = 0.01; end
    if ~isfield(Params,'gThresh'), Params.gThresh = 0.64; end
    if ~isfield(Params,'rThresh'), Params.rThresh = 10; end
    if ~isfield(Params,'wr'), Params.wr = 1; end
    if ~isfield(Params,'ws'), Params.ws = 1; end
    if ~isfield(Params,'wd'), Params.wd = 5; end
    if ~isfield(Params,'wn'), Params.wn = 0; end
    if ~isfield(Params,'smoothing'), Params.smoothing = true; end
    if ~isfield(Params,'gaussSig'), Params.gaussSig = 4; end%4/(norm(rSat)/35); end
end

lthresh = Params.lthresh;
dThreshUp = Params.dThreshUp;
dThreshLo = Params.dThreshLo;
ceSml = Params.ceSml;
ceLrg = Params.ceLrg;
sThresh = Params.sThresh;
gThresh = Params.gThresh;
rThresh = Params.rThresh;
wr = Params.wr;
wd = Params.wd;
wn = Params.wn;
ws = Params.ws;
smoothing = Params.smoothing;
gaussSig = Params.gaussSig;
szMin = 12;
    
% if nanmean(nanmean(Im)) < 0.02, Ldxn = -Ldxn; end
dat = transPt32(-rAst',rSat,R,camera);
Ldxn = -dat./norm(dat);
a = acosd(dot(rSat,-rAst)/(norm(rSat)*norm(-rAst)));
if a > 90
    Ldxn = -Ldxn;
end
[nR,nC] = size(Im);
Lnrm = Ldxn./norm(Ldxn);

Sx = Lnrm(1); % Sx and Sy need to point TOWARDS the light
Sy = Lnrm(2); 

% Find corner to use for intersect of light line
if Ldxn(1) >= 0 && Ldxn(2) >= 0
    % Bottom Left Corner
    int = [1,1];
elseif Ldxn(1) >= 0 && Ldxn(2) < 0
    % Top Left Corner
    int = [nR,1];
elseif Ldxn(1) < 0 && Ldxn(2) < 0
    % Top Right Corner
    int = [nR,nC];
else % Ldxn(1) < 0 && Ldxn(2) >= 0
    % Bottom Right Corner
    int = [1,nC];
end
m = -1/(Lnrm(2)/Lnrm(1));
A = [-m 1 -(int(1)-m*int(2))]; % (a,b,c) of a*col+b*(nR-row+1)+c=0

%% Edge Detection
if smoothing
    % ImBlurred = imgaussfilt(Im,gaussSig);
    ImCont = imadjust(Im);
    kernel = ones(gaussSig) / gaussSig^2;
    ImBlurred = imfilter(ImCont, kernel);
else
    ImBlurred = Im;
end
[Gx,Gy] = imgradientxy(ImBlurred);
Esml = edge(ImBlurred,'canny',[0 ceSml]); % Small, sharp edges
Elrg = edge(ImBlurred,'canny',[ceSml ceLrg]); % Small, sharp edges
Eog = Esml;
Eog(Elrg) = 1;
% Remove Edges w/ Incorrect Gradient
Gck = ((Gx.*Sx + Gy.*Sy)./(Gx.^2 + Gy.^2).^(1/2)) > gThresh;
Egr = Eog.*Gck;

%% Identify Lit and Shaded Edges
CC = bwconncomp(Egr);
Ecc = Egr;

% Remove small edges
for o = 1:CC.NumObjects
    if length(CC.PixelIdxList{o}) < szMin
        Ecc(CC.PixelIdxList{o}) = 0;
    end
end

CC = bwconncomp(Ecc);
    
for o = 1:CC.NumObjects
    CC.PixelSubList{o} = zeros(length(CC.PixelIdxList{o}),2);
    for p = 1:length(CC.PixelIdxList{o})
        [r,c] = ind2sub([nR,nC],CC.PixelIdxList{o}(p));
        CC.PixelSubList{o}(p,:) = [r,c];
    end
end

if CC.NumObjects == 0
    craters = [];
    CC = [];
    ePar = [];
    crit = [];
    return
end

CC.DistToLight = CC.PixelIdxList;
CC.CenterPx = zeros(1,CC.NumObjects);
CC.ClosestPx = zeros(1,CC.NumObjects);
CC.FarthestPx = zeros(1,CC.NumObjects);

% Split into Shaded and Lit Sides
LitSides = zeros(1,CC.NumObjects);
Ess = Egr;
Els = Egr;
rmv = [];
for o = 1:CC.NumObjects
    px = CC.PixelIdxList{o};
    for p = 1:length(px)
        r = CC.PixelSubList{o}(p,1);
        c = CC.PixelSubList{o}(p,2);
        CC.DistToLight{o}(p) = abs(A(1)*c + A(2)*(nR-r+1) + A(3))/...
            sqrt(A(1)^2 + A(2)^2);
    end
    [dc,pc] = min(CC.DistToLight{o});
    [df,pf] = max(CC.DistToLight{o});
    pm = round(length(px)/2);
    CC.CenterPx(o) = pm;
    CC.ClosestPx(o) = pc;
    CC.FarthestPx(o) = pf;
    r = CC.PixelSubList{o}(pm,1);
    c = CC.PixelSubList{o}(pm,2);
    Pxm = [c;nR-r+1];
    r = CC.PixelSubList{o}(1,1);
    c = CC.PixelSubList{o}(1,2);
    Px1 = [c;nR-r+1];
    r = CC.PixelSubList{o}(end,1);
    c = CC.PixelSubList{o}(end,2);
    Px2 = [c;nR-r+1];
    dm = CC.DistToLight{o}(pm);
    de1 = CC.DistToLight{o}(1);
    de2 = CC.DistToLight{o}(end);
    Gm = [Gx(px(pm));Gy(px(pm))];
    G1 = [Gx(px(1));Gy(px(1))];
    G2 = [Gx(px(end));Gy(px(end))];
    tmp = dot((Pxm-Px1),G1)/(norm(Pxm-Px1)*norm(G1)) + ...
          dot((Pxm-Px2),G2)/(norm(Pxm-Px2)*norm(G2));
    % sort
    if tmp <= 0
        % Shaded Side
        Els(CC.PixelIdxList{o}) = 0;
    else
        % Lit Side
        LitSides(o) = 1;
        Ess(CC.PixelIdxList{o}) = 0;
    end
end

CC.LitSides = LitSides;

% For speed...
ImageSize = CC.ImageSize;
NumObjects = CC.NumObjects;
PixelIdxList = CC.PixelIdxList;
PixelSubList = CC.PixelSubList;
DistToLight = CC.DistToLight;
CenterPx = CC.CenterPx;
ClosestPx = CC.ClosestPx;
FarthestPx = CC.FarthestPx;
LitSides = CC.LitSides;

%% Pair off Edges
mtchd = zeros(NumObjects,1);
ePar = [(1:NumObjects)' zeros(NumObjects,6)];
crit = [0,0,0,0]; % l_rto, dist, cone, slope
for es = 1:(NumObjects-1)
    % test: 843
    ls = length(PixelIdxList{es});
    [rms,cms] = ind2sub([nR,nC],PixelIdxList{es}(round(ls/2)));
    [r1s,c1s] = ind2sub([nR,nC],PixelIdxList{es}(1));
    [r2s,c2s] = ind2sub([nR,nC],PixelIdxList{es}(end));
    best = [0 100000 0 100000 0 100000];
    
    % Cone of Light
    Pcen = PixelSubList{es}(CenterPx(es),[2,1]);
    if LitSides(es)
        Pext = Pcen+Ldxn;
    else
        Pext = Pcen-Ldxn;
    end
    [mN,bN] = ptSlope(Pcen,Pcen+Ldxn); % Normal
    [mT,bT] = ptSlope(Pcen,Pcen+fliplr(Ldxn)); % Tangent
    Pcw = (Pcen' + rot2Dd(15)*Ldxn')';
    Pccw = (Pcen' + rot2Dd(-15)*Ldxn')';
    [mCW,bCW] = ptSlope(Pcen,Pcw);
    [mCCW,bCCW] = ptSlope(Pcen,Pccw);
    sns = sign([Pext(2)-(mT*Pext(1)+bT), Pext(2)-(mCW*Pext(1)+bCW), Pext(2)-(mCCW*Pext(1)+bCCW)]);
    
    % Plot Cone:
    %{
        figure; imshow(cn)
    %}
    
    % Random match order
    el_list = (es+1):NumObjects;
    el_list = el_list(randperm(length(el_list)));
    
    for el = el_list
        % test: 871
        mtch = true;
        % any(PixelSubList{es}(:,2) == 1404) && any(PixelSubList{es}(:,1) == 1026) && any(PixelSubList{el}(:,2) == 1454) && any(PixelSubList{el}(:,1) == 982)
        PSL = PixelSubList{el};
        PIL = PixelIdxList{el};
        % Plot edges:
        %{
                figure(); imshow(Im); hold on
                plot(PixelSubList{el}(:,2),PixelSubList{el}(:,1),'y')
                plot(PixelSubList{es}(:,2),PixelSubList{es}(:,1),'g')
        %}
        % Length Ratio Check
        ll = size(PSL,1);
        rto = ll/ls;
        if 1/lthresh > rto || rto > lthresh
            mtch = false;
            crit(1) = crit(1) + 1;
        end
        % Distance Between Centers
        if mtch
            rel = round(ll/2);
            ml = PSL(rel,:);
            rml = ml(1);
            cml = ml(2);
            d = sqrt((cms-cml)^2+(rms-rml)^2);
            if (ll >= ls) && (d > dThreshUp*ll || d < dThreshLo*ls)
                mtch = false;
                crit(2) = crit(2) + 1;
            elseif (ll < ls) && (d > dThreshUp*ls || d < dThreshLo*ll)
                mtch = false;
                crit(2) = crit(2) + 1;
            end
            %                 if abs(d) < 2*szMin
            %                     mtch = false;
            %                 end
        end
        % Within 30 degree cone
        if mtch
            cnCK = zeros(1,ll);
            for idx = 1:ll
                r = PSL(idx,1);
                c = PSL(idx,2);
                snsCK = sign([r-(mT*c+bT), r-(mCW*c+bCW), r-(mCCW*c+bCCW)]);
                if all(snsCK == sns)
                    cnCK(idx) = 1;
                end
            end
            if sum(cnCK)/ll < 0.6
                mtch = false;
                crit(3) = crit(3) + 1;
            end
            %                 PcenC = PSL(CenterPx(el),[2,1]);
            %                 snsL = sign([PcenC(2)-(mT*PcenC(1)+bT), PcenC(2)-(mCW*PcenC(1)+bCW), PcenC(2)-(mCCW*PcenC(1)+bCCW)]);
            %                 Pcen1 = PSL(1,[2,1]);
            %                 snsC = sign([Pcen1(2)-(mT*Pcen1(1)+bT), Pcen1(2)-(mCW*Pcen1(1)+bCW), Pcen1(2)-(mCCW*Pcen1(1)+bCCW)]);
            %                 Pcen2 = PSL(end,[2,1]);
            %                 snsF = sign([Pcen2(2)-(mT*Pcen2(1)+bT), Pcen2(2)-(mCW*Pcen2(1)+bCW), Pcen2(2)-(mCCW*Pcen2(1)+bCCW)]);
            %                 if any(snsC~=sns) || ~( all(snsL==sns) || all(snsC==sns) || all(snsF==sns))
            %                     mtch = false;
            %                 end
        end
        % Linearly Increasing Intensity
        if mtch
            x = [cml cms];
            y = [rml rms];
            c = improfile(Im,x,y,10);
            [s,res] = polyfit(1:length(c),c',1);
            if LitSides(es)
                if (s(1) < 0) || (abs(s(1)) < sThresh) %|| (res.normr > rThresh)
                    mtch = false;
                    crit(4) = crit(4) + 1;
                end
            else
                if (s(1) > 0) || (abs(s(1)) < sThresh) %|| (res.normr > rThresh)
                    mtch = false;
                    crit(4) = crit(4) + 1;
                end
            end
            cst = wr*abs(rto-1) + ws/abs(s(1)) + wd*d/ll + wn*res.normr;
        end
        % Match Analysis
        if mtch && (cst < best(6))
            best = [el rto s(1) res.normr d cst];
        end
    end
    ePar(find(ePar(:,1)==es),2:7) = best;
%     ePar(find(ePar(:,1)==el),:) = [];
end

%% Remove Double Pairings
ePar(ePar(:,2)==0,:) = [];
for idx = unique(ePar(:,2)')
    ct = sum(ePar(:,2) == idx);
    if ct > 1
        k = find(ePar(:,2) == idx);
        cst = ePar(k,7);
        [~,bst] = min(cst);
        k(bst) = [];
        ePar(k,:) = [];
    end
end

%% Ellipse Fitting
nCraters = size(ePar,1);
craters = zeros(nCraters,6);

% Direct Method A. W. Fitzgibbon, M. Pilu, R. B. Fishe
% https://www.mathworks.com/matlabcentral/fileexchange/22684-ellipse-fit-direct-method
for cr = 1:nCraters
    YX = [CC.PixelSubList{ePar(cr,1)};CC.PixelSubList{ePar(cr,2)}];
    x = YX(:,2);
    y = YX(:,1);
    if size(YX,1) > 5
        [~,B] = EllipseDirectFit([YX(:,2),YX(:,1)]);        
        craters(cr,1:6) = B; % [x0 y0 a b alpha score]
    end
end

%% Eliminate Overlapping Craters
elim = zeros(nCraters,1);
for c1 = 1:(nCraters-1)
    % Find overlapping craters
    a = 1.2*craters(c1,3);
    lims = [craters(c1,1)-a,craters(c1,1)+a,craters(c1,2)-a,craters(c1,2)+a];
    ovlp = (craters((c1+1):end,1)>lims(1)).*(craters((c1+1):end,1)<lims(2)).*...
           (craters((c1+1):end,2)>lims(3)).*(craters((c1+1):end,2)<lims(4));
    if sum(ovlp) > 0
        % Keep the smallest crater
        c2s = find(ovlp)+c1;
        kp = c1;
        elim(c1) = 1;
        for c2 = c2s
            elim(c2) = 1;
            if craters(c2,3) < craters(kp,3)
                kp = c2;
            end
        end
        elim(kp) = 0;
    end
end

craters(find(elim),:) = [];

end

function [m,b] = ptSlope(pt1,pt2)
% pt = [x,y]
m = (pt1(2)-pt2(2))/(pt1(1)-pt2(1));
b = pt1(2) - m*pt1(1);
end