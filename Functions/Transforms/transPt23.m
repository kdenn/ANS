function ftACI = transPt23(points2D,rSat,R,f,vACI,camera)
% Convert positions from the 2D image pixel frame to 3D ACI ecliptic
% INPUT:
    % points2D: pixel col/row coordinates [c r] [nx2]
    % rSat: position of satellite in ACI Ecliptic [3x1]
    % R: 3x3 rotation matrix from ACI Ecliptic to Camera (see rotWldCam)
    % f: faces defined by vertex numbers [nfx3] [v1 v2 v3]
    % vACI: matrix of vertices [nvx3] [x y z] in ACI Ecliptic Plane
    % camParams: array of camera parameters
        % [foc, pxNum, pxSze, radDist, imPxAct]
% OUTPUT:
    % ftACI: [x y z] feature positions in ACI frame [nx3]
    
np = size(points2D,1);
nv = size(vACI,1);
nf = size(f,1);
X = zeros(np,3);
ftACI = zeros(np,3);
dv = norm(vACI(f(1,1),:)-vACI(f(1,2),:));

%% Get Camera Intrinsics
% AA273 AR lec 5
pxNum = camera.pxNum; % pixels in CCD pixel array cxr
t = R*(-rSat);
P = camera.A*[R t];
Pi = P'*(P*P')^(-1);

%% Camera Transformation
% Convert px coordinates back to non-square px coordinates
points2D(:,1) = pxNum(1) - points2D(:,1);
% Get Second Point of Ray from Camera to Feature
for i = 1:np
    ph = [points2D(i,[1,2])'; 1];
    PWh = Pi*ph;
    X(i,:) = (PWh(1:3)./PWh(4))';
end

%% Calculate Body Intersection
db = zeros(nv,1);
C = rSat';
for i = 1:np
    Dxn = X(i,:)-C;
    t_best = inf;
    for j = 1:nv
        d = abs(cross(vACI(j,:)-X(i,:),C))/abs(C);
        % http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        if d <= 3*dv
            db(j) = 1;
        end
    end
    for j = 1:nf
        if all(db(f(j,:)))
            [flag,u,v,t] = rayTriangleIntersection(C,Dxn,vACI(f(j,1),:),vACI(f(j,2),:),vACI(f(j,3),:));
            if flag && t <= t_best
                ftACI(i,1:3) = C + t*Dxn;
                t_best = t;
            end
        end
    end
end

end

function [flag, u, v, t] = rayTriangleIntersection (o, d, p0, p1, p2)
% Ray/triangle intersection using the algorithm proposed by Möller and Trumbore (1997).
%
% Input:
%    o : origin.
%    d : direction.
%    p0, p1, p2: vertices of the triangle.
% Output:
%    flag: (0) Reject, (1) Intersect.
%    u,v: barycentric coordinates.
%    t: distance from the ray origin.
% Author: 
%    Jesus Mena

    epsilon = 0.00001;

    e1 = p1-p0;
    e2 = p2-p0;
    q  = cross(d,e2);
    a  = dot(e1,q); % determinant of the matrix M

    if (a>-epsilon && a<epsilon) 
        % the vector is parallel to the plane (the intersection is at infinity)
        [flag, u, v, t] = deal(0,0,0,0);
        return;
    end;
    
    f = 1/a;
    s = o-p0;
    u = f*dot(s,q);
    
    if (u<0.0)
        % the intersection is outside of the triangle
        [flag, u, v, t] = deal(0,0,0,0);
        return;          
    end;
    
    r = cross(s,e1);
    v = f*dot(d,r);
    
    if (v<0.0 || u+v>1.0)
        % the intersection is outside of the triangle
        [flag, u, v, t] = deal(0,0,0,0);
        return;
    end;
    
    t = f*dot(e2,r); % verified! 
    flag = 1;
    return;
end