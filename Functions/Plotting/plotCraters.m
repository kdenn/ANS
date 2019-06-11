function plotCraters(craters,var1,var2,var3)

if nargin == 1
    % plotCraters(crater)
    for c = 1:size(craters,1)
        ellipse(craters(c,3),craters(c,4),craters(c,5),craters(c,1),craters(c,2),'r')
    end
elseif nargin == 2
    % plotCraters(craters,color)
    for c = 1:size(craters,1)
        ellipse(craters(c,3),craters(c,4),craters(c,5),craters(c,1),craters(c,2),var1)
    end
elseif nargin == 3
    % plotCraters(craters,CC,ePar)
    for c = 1:size(craters,1)
        plot(var1.PixelSubList{var2(c,1)}(:,2),var1.PixelSubList{var2(c,1)}(:,1),'y')
        plot(var1.PixelSubList{var2(c,2)}(:,2),var1.PixelSubList{var2(c,2)}(:,1),'c')
        ellipse(craters(c,3),craters(c,4),craters(c,5),craters(c,1),craters(c,2))
    end
else 
    % plotCraters(craters,CC,ePar,color)
    for c = 1:size(craters,1)
        plot(var1.PixelSubList{var2(c,1)}(:,2),var1.PixelSubList{var2(c,1)}(:,1),'y')
        plot(var1.PixelSubList{var2(c,2)}(:,2),var1.PixelSubList{var2(c,2)}(:,1),'c')
        ellipse(craters(c,3),craters(c,4),craters(c,5),craters(c,1),craters(c,2),var3)
    end
end
end