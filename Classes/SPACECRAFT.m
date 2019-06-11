classdef SPACECRAFT
    % Spacecraft with a camera orbiting an asteroid
    % Kaitlin Dennison - SLAB - Spring 2019
    properties
        name
        camera                  % camera object
    end
    properties (SetAccess = protected)
        tEpoch                  % Julian date of epoch (days)
        oeACI                   % ACI orbital elements [a e i Om w M_0]' (km - rad)
        mu                      % gravitational parameter of the CB
    end
    methods
        function sc = SPACECRAFT(JD_0,oeACI,mu,name)
            % Initialize the spacecraft
            % oeACI = [a e i Om w M_0]'
            sc.tEpoch = JD_0;
            sc.oeACI = oeACI;
            sc.mu = mu;
            if nargin > 3
                sc.name = name;
            end
        end
        function [rACI,vACI] = getRV(sc,JD)
            % Get the ACI position and velocity of the SC at JD
            M = mod(sc.M_0 + sc.n*(JD-sc.JD_epoch)*60*60*24,2*pi);
            Ec = iterMtoEc(M,sc.e,1E-8);
            nu = acos((cos(Ec)-sc.e)/(1-sc.e* cos(Ec)));
            if Ec > pi
                nu = 2*pi-nu;
            end
            [rACI,vACI] = OEtoRV(sc.a,sc.e,sc.i,sc.Om,sc.w,nu,sc.mu);
        end
    end
end