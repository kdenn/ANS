function plotMatch(Im,craters,ftDB2D,ftDB,corr,CI)
% Plot the image and how the database was correlated
% INPUT:
%   Im - the image
%   craters - output from detectCraters
%   ftDB2D - 2D projections of database
%       ftDB2D = transPt32(ftACI,rSatACI,rotACItoCam,camParams)
%   ftDB - the feature database with current 2D cov
%   corr - output from corrFeat
%   CI - confidence interval for mah
%       https://people.richland.edu/james/lecture/m170/tbl-chi.html
%--------------------------------------------------------------------------

figure(); hold on
imshow(Im); hold on
% scatter(craters(:,1),craters(:,2),'^c') % New Craters
scatter(ftDB2D(corr(:,1),1),ftDB2D(corr(:,1),2),'*g') % Correlated database craters
for c = 1:size(ftDB2D,1)
    % plot the 2D cov ellipses
    % http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
    if ~any(corr(:,1)==c)
        continue
    end
    [V,D] = eig(getDBfeat(ftDB,c,'2Dcov'));
    [a,i] = max(diag(D));
    a = 2*sqrt(CI*a);
    if i == 1
        b = 2*sqrt(CI*D(2,2));
    else
        b = 2*sqrt(CI*D(1,1));
    end
    alpha = atan(V(2,i)/V(1,i));
    ellipse(a,b,alpha,ftDB2D(c,1),ftDB2D(c,2),'g')
end
scatter(corr(:,2),corr(:,3),'+r') % Correlated new craters
for c = 1:size(craters,1)
    ellipse(craters(c,3),craters(c,4),craters(c,5),craters(c,1),craters(c,2))
end
hold off
hold off

end