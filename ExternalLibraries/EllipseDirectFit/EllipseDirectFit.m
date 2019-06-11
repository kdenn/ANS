function [A,B] = EllipseDirectFit(XY)
%
%  Direct ellipse fit, proposed in article
%    A. W. Fitzgibbon, M. Pilu, R. B. Fisher
%     "Direct Least Squares Fitting of Ellipses"
%     IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
%
%  Our code is based on a numerically stable version
%  of this fit published by R. Halir and J. Flusser
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: A = [a b c d e f]' is the vector of algebraic 
%             parameters of the fitting ellipse:
%             ax^2 + bxy + cy^2 +dx + ey + f = 0
%             the vector A is normed, so that ||A||=1
%
%  This is a fast non-iterative ellipse fit.
%
%  It returns ellipses only, even if points are
%  better approximated by a hyperbola.
%  It is somewhat biased toward smaller ellipses.
%
centroid = mean(XY);   % the centroid of the data set

D1 = [(XY(:,1)-centroid(1)).^2, (XY(:,1)-centroid(1)).*(XY(:,2)-centroid(2)),...
      (XY(:,2)-centroid(2)).^2];
D2 = [XY(:,1)-centroid(1), XY(:,2)-centroid(2), ones(size(XY,1),1)];
S1 = D1'*D1;
S2 = D1'*D2;
S3 = D2'*D2;
if det(S3) == 0
    T = -pinv(S3)*S2';
else
    T = -inv(S3)*S2';
end
M = S1 + S2*T;
M = [M(3,:)./2; -M(2,:); M(1,:)./2];
[evec,eval] = eig(M);
cond = 4*evec(1,:).*evec(3,:)-evec(2,:).^2;
A1 = evec(:,find(cond>0));
A = [A1; T*A1];
A4 = A(4)-2*A(1)*centroid(1)-A(2)*centroid(2);
A5 = A(5)-2*A(3)*centroid(2)-A(2)*centroid(1);
A6 = A(6)+A(1)*centroid(1)^2+A(3)*centroid(2)^2+...
     A(2)*centroid(1)*centroid(2)-A(4)*centroid(1)-A(5)*centroid(2);
A(4) = A4;  A(5) = A5;  A(6) = A6;
A = A/norm(A);

% Get Fit
x = XY(:,1);
y = XY(:,2);
D = [x.*x x.*y y.*y x y ones(size(x))];
score = std(A'*D'); % Standard deviation of the ellipse fit

% ax^2 + bxy + cy^2 +dx + ey + f = 0

% Get the ellipse parameters
a = A(1); b = A(2)/2; c = A(3); d = A(4)/2; f = A(5)/2; g = A(6);

center(1) = (c*d - b*f)/(b^2-a*c);
center(2) = (a*f - b*d)/(b^2-a*c);

sem(1) = sqrt( 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g) / ((b^2-a*c)*(sqrt((a-c)^2+4*b^2)-(a+c))));
sem(2) = sqrt( 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g) / ((b^2-a*c)*(-sqrt((a-c)^2+4*b^2)-(a+c))));

if b == 0 && a < c
phi = 0;
elseif b == 0 && a > c
phi = 0.5*pi;
elseif b ~= 0 && a < c
phi = 0.5* acot((a-c)/(2*b));
else
phi = 0.5*pi + 0.5* acot((a-c)/(2*b));
end

B = [center,sem,phi,score];

% Plot Fit:
%{
    hold on
    %Convert the A to str
    a = num2str(A(1));
    b = num2str(A(2));
    c = num2str(A(3));
    d = num2str(A(4));
    e = num2str(A(5));
    f = num2str(A(6));
    %Points
    plot(XY(:,1),XY(:,2),'rx')
    %Equation
    eqt= ['(',a, ')*x^2 + (',b,')*x*y + (',c,')*y^2 + (',d,')*x+ (',e,')*y + (',f,')'];
    xmin=0.7*min(XY(:,1));
    xmax=1.3*max(XY(:,2));
    ezplot(eqt,[xmin,xmax])
    scatter(XY(:,1),XY(:,2))
    %Ellipse Paramers
    ellipse(sem(1),sem(2),phi,center(1),center(2))
    hold off
%}

end  %  EllipseDirectFit


