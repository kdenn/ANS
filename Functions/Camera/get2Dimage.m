function Im = get2Dimage(f,v,rAst,rSat,R,camera)
% Get a 2D image from the vertex data
% INPUT
    % f: faces defined by vertex numbers [nfx3] [v1 v2 v3]
    % v: matrix of vertices [nvx3] [x y z] in ACI Equatorial Plane
    % rAst: position of satellite in SCI Equatorial [3x1]
    % rSat: position of satellite in ACI Equatorial [3x1]
    % R: 3x3 rotation matrix from ACI Equatorial to Camera (see rotWldCam)
    % camera: camera object describing all the intrinsics
% OUTPUT
    % Im: the 2D image of the asteroid from the sat/cam perspective

nv = size(v,1);
nf = size(f,1);
v_cam = zeros(nv,2);
rSatAst = -rSat;
[FaceVertexCData,fOrder] = getLighting(f,v,rAst,rSat);
FaceVertexCData = 0.8.*FaceVertexCData;

%% Camera Constants
% https://nssdc.gsfc.nasa.gov/nmc/experimentDisplay.do?id=1996-008A-01
pxNum = camera.pxNum; % pixels in CCD pixel array rxc
t = R*(-rSat);
K = camera.A*[R t];

%% Convert from World to 2D Camera
% AA273 AR lec 5 s13
for i = 1:nv
    Pw = v(i,:)';
    Pwh = [Pw; 1];
    ph = K*Pwh;
    p = ph(1:2)./ph(3);
    v_cam(i,:) = p';
end

%% Determine which faces are actually visible
sc = 3;
vizFaces = getViz(f,v,rSatAst);
xs = ([v_cam(f(:,1),1) v_cam(f(:,2),1) v_cam(f(:,3),1)]).*sc;
ys = ([v_cam(f(:,1),2) v_cam(f(:,2),2) v_cam(f(:,3),2)]).*sc;

%% Create the Image 
pxNumSc = pxNum.*sc;
Im = zeros(pxNum);
for fi = fOrder'
    if vizFaces(fi)
        xr = floor(min(xs(fi,:))):ceil(max(xs(fi,:)));
        yr = floor(min(ys(fi,:))):ceil(max(ys(fi,:)));
        xr(xr>pxNumSc(1)) = [];
        xr(xr<1) = [];
        yr(yr>pxNumSc(2)) = [];
        yr(yr<1) = [];
        if ~isempty(xr) && ~isempty(yr)
            [Xq,Yq] = meshgrid(xr,yr);
            [in,on] = inpolygon(Xq,Yq,xs(fi,:),ys(fi,:));
            if any(any(in))
                Im(Xq(in),Yq(in)) = FaceVertexCData(fi,1);
            end
            if any(any(on))
                Im(Xq(on),Yq(on)) = FaceVertexCData(fi,1);
            end
        end
    end
end

Im = imresize(Im,pxNum);
Im = rot90(Im);
Im = flip(Im,1);
Im = flip(Im,2);

if camera.radDist ~= 0, Im = lensdistort(Im,camera.radDist); end

end

function vizFaces = getViz(f,v,rSatAst)
nf = size(f,1);
vizFaces = zeros(nf,1);
theta = zeros(nf,1);
normals = zeros(nf,3);
centers = zeros(nf,3);
rFacSun = zeros(nf,3);
if size(rSatAst,1) == 3
    rSatAst = rSatAst';
end

for i = 1:nf
    p1 = v(f(i,1),:);
    p2 = v(f(i,2),:);
    p3 = v(f(i,3),:);
    U = p2 - p1;
    V = p3 - p1;
    N = cross(U,V);
    normals(i,:) = N;
    centers(i,:) = [v(f(i,1),1)+v(f(i,2),1)+v(f(i,3),1), ...
                    v(f(i,1),2)+v(f(i,2),2)+v(f(i,3),2), ...
                    v(f(i,1),3)+v(f(i,2),3)+v(f(i,3),3)]./3;
    rFacSun(i,:) = -(rSatAst+centers(i,:));
    theta(i) = acos(dot(normals(i,:)./norm(normals(i,:)),rFacSun(i,:)./norm(rFacSun(i,:))));
    if theta(i) > pi/2
        vizFaces(i) = 0;
    else
        vizFaces(i) = 1;
    end
end
vizFaces = logical(vizFaces);
end

function [FaceVertexCData, fOrder] = getLighting(f,v,rSunAst,rSat)
% Get the light intesity of each face
% INPUT
%   f: faces defined by vertex numbers [nfx3] [v1 v2 v3]
%   v: matrix of vertices [nvx3] [x y z] in ACI Plane
% normals: https://www.khronos.org/opengl/wiki/Calculating_a_Surface_Normal
% centers: http://mathforum.org/library/drmath/view/54899.html
nf = size(f,1);
rho = 1;
ISun = 1;
Iface = zeros(nf,1);
theta = zeros(nf,1);
normals = zeros(nf,3);
rAstFac = zeros(nf,3);
rFacSun = zeros(nf,3);
rSatFac = zeros(nf,3);
FaceVertexCData = ones(nf,3);
if size(rSunAst,1) == 3
    rSunAst = rSunAst';
end

for i = 1:nf
    p1 = v(f(i,1),:);
    p2 = v(f(i,2),:);
    p3 = v(f(i,3),:);
    U = p2 - p1;
    V = p3 - p1;
    N = cross(U,V);
    normals(i,:) = N;
    rAstFac(i,:) = [v(f(i,1),1)+v(f(i,2),1)+v(f(i,3),1), ...
                    v(f(i,1),2)+v(f(i,2),2)+v(f(i,3),2), ...
                    v(f(i,1),3)+v(f(i,2),3)+v(f(i,3),3)]./3; %ACI
    rFacSun(i,:) = -(rSunAst+rAstFac(i,:));
    theta(i) = acos(dot(normals(i,:)./norm(normals(i,:)),rFacSun(i,:)./norm(rFacSun(i,:))));
    if theta(i) > pi/2 
        Iface(i) = 0;
    else
        Iface(i) = rho*ISun*cos(theta(i));
    end
    FaceVertexCData(i,:) = ones(1,3).*Iface(i);
    rSatFac(i,:) = (rSat'+rAstFac(i,:));
end

fDist = (rSatFac(:,1).^2+rSatFac(:,2).^2+rSatFac(:,3).^2).^(1/2);
[~,fOrder] = sort(fDist,'ascend');

end


















