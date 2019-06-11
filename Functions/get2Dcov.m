function [M,xp] = get2Dcov(xACAF,P,hmeas)
% From "Probabilistic Robotics" by Thrun p.66

xACAF = xACAF';

% Parameters
L = 2;
n = length(xACAF);
k = 0;
b = 2;
a = sqrt((L+n)/(n+k));

% State Sigma Points
nP = 2*n+1;
C = sqrt(n+L);
A = C*chol(P)';
Xi = xACAF(:,ones(1,n));
X = [xACAF Xi+A Xi-A]; 

% Weighting Matrices
Wm=[L/(n+L) 0.5/(n+L)+zeros(1,2*n)];
Wc=Wm;
Wc(1)=Wc(1)+(1-a^2+b);

% Measurement Sigma Points
Z = zeros(2,nP);
for i = 1:nP
    Z(:,i) = hmeas(X(:,i));
end

% Map Points
xp = zeros(2,1);
for i = 1:nP
    xp = xp + Wm(:,i)*Z(:,i);
end
M = zeros(2);
for i = 1:nP
    M = M + Wc(:,1)*(Z(:,i)-xp)*(Z(:,i)-xp)';
end

end

%{
function M = get2Dcov(xACI,R,Pvec,rSat,camParams)
rx = xACI(1);
ry = xACI(2);
rz = xACI(3);
P_ACI = diag(R*Pvec);

%% Camera Constants
foc = camParams(1); % [mm] focal length
pxNum = camParams(2:3); % pixels in CCD pixel array rxc
pxSze = camParams(4:5); % [mm] wxh of one px
imSze = pxSze.*pxNum([2,1]); % [mm] size of CCD array
k = 1./pxSze; % [px/mm]

%% Get Camera Intrinsics
% AA273 AR lec 5
x0 = imSze(1)/2;
y0 = imSze(2)/2;
t = R*(-rSat);
fx = foc*k(1);
fy = foc*k(2);
kx = k(1)*x0;
ky = k(2)*y0;
tx = t(1);
ty = t(2);
tz = t(3);

%% Compute the Jacobian
J = zeros(2,3);
J(1,1) = (R(1,1)*fx + R(3,1)*kx)/(tz + R(3,1)*rx + R(3,2)*ry + R(3,3)*rz) - (R(3,1)*(fx*tx + kx*tz + rx*(R(1,1)*fx + R(3,1)*kx) + ry*(R(1,2)*fx + R(3,2)*kx) + rz*(R(1,3)*fx + R(3,3)*kx)))/(tz + R(3,1)*rx + R(3,2)*ry + R(3,3)*rz)^2;
J(1,2) = (R(1,2)*fx + R(3,2)*kx)/(tz + R(3,1)*rx + R(3,2)*ry + R(3,3)*rz) - (R(3,2)*(fx*tx + kx*tz + rx*(R(1,1)*fx + R(3,1)*kx) + ry*(R(1,2)*fx + R(3,2)*kx) + rz*(R(1,3)*fx + R(3,3)*kx)))/(tz + R(3,1)*rx + R(3,2)*ry + R(3,3)*rz)^2;
J(1,3) = (R(1,3)*fx + R(3,3)*kx)/(tz + R(3,1)*rx + R(3,2)*ry + R(3,3)*rz) - (R(3,3)*(fx*tx + kx*tz + rx*(R(1,1)*fx + R(3,1)*kx) + ry*(R(1,2)*fx + R(3,2)*kx) + rz*(R(1,3)*fx + R(3,3)*kx)))/(tz + R(3,1)*rx + R(3,2)*ry + R(3,3)*rz)^2;
J(2,1) = (R(2,1)*fy + R(3,1)*ky)/(tz + R(3,1)*rx + R(3,2)*ry + R(3,3)*rz) - (R(3,1)*(fy*ty + ky*tz + rx*(R(2,1)*fy + R(3,1)*ky) + ry*(R(2,2)*fy + R(3,2)*ky) + rz*(R(2,3)*fy + R(3,3)*ky)))/(tz + R(3,1)*rx + R(3,2)*ry + R(3,3)*rz)^2;
J(2,2) = (R(2,2)*fy + R(3,2)*ky)/(tz + R(3,1)*rx + R(3,2)*ry + R(3,3)*rz) - (R(3,2)*(fy*ty + ky*tz + rx*(R(2,1)*fy + R(3,1)*ky) + ry*(R(2,2)*fy + R(3,2)*ky) + rz*(R(2,3)*fy + R(3,3)*ky)))/(tz + R(3,1)*rx + R(3,2)*ry + R(3,3)*rz)^2;
J(2,3) = (R(2,3)*fy + R(3,3)*ky)/(tz + R(3,1)*rx + R(3,2)*ry + R(3,3)*rz) - (R(3,3)*(fy*ty + ky*tz + rx*(R(2,1)*fy + R(3,1)*ky) + ry*(R(2,2)*fy + R(3,2)*ky) + rz*(R(2,3)*fy + R(3,3)*ky)))/(tz + R(3,1)*rx + R(3,2)*ry + R(3,3)*rz)^2;
M = abs(J*P_ACI*J');
end
%}