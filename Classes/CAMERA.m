classdef CAMERA
    % Camera intrinsics and image propagation
    % Kaitlin Dennison - SLAB - Spring 2019
    properties
        name
    end
    properties (SetAccess = protected)
        f                       % focal length (mm)
        pxNum                   % sensor size [wxh] (px)
        pxSze                   % size of a pixel [wxh] (micron)
        A                       % instrinsic matrix [3x3]
        radDist = 0             % radial distortion
    end
    methods
        function c = CAMERA(f,pxNum,pxSze,name,radDist)
            % Camera constructor
            c.f = f;
            c.pxNum = pxNum;
            c.pxSze = pxSze;
            k = 1000./pxSze; % (px/mm)
            imSze = pxSze.*pxNum.*0.001; % [wxh] (mm)
            x0 = imSze(1)/2;
            y0 = imSze(2)/2;
            c.A = [f*k(1) 0      k(1)*x0;
                   0      f*k(2) k(2)*y0;
                   0      0      1];
            if nargin > 3
                c.name = name;
            end
            if nargin == 5
                c.radDist = radDist;
            end
        end
    end
end