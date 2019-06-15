classdef ASTEROID
    % Short description (fit in one line)
    %{
    -----------------------------------------------------------------------
    DESCRIPTION:
        Longer description here in paragraph form. All text should fit
        within the standard MATLAB width unless the line includes a string.
    -----------------------------------------------------------------------
    REFERENCES:
        - Citations here, preferrably IEEE format
    -----------------------------------------------------------------------
    NOTES:
        - Text text text text text text text text text text text text text
          text text text text text text text text text text text text
        - Text text text text text text text text text text text text text
          text text text text text text text text text text text text
    -----------------------------------------------------------------------
    AUTHOR: Kaitlin Dennison
    -----------------------------------------------------------------------
    COPYTRIGHT: 2019 SLAB Group
    -----------------------------------------------------------------------
    TIMESTAMPS:
        - 14-Jun-2019: creation (KD)
    -----------------------------------------------------------------------
    %}
    
    properties (SetAccess = private)
        name = 'Eros'               % str, name of the asteroid
        filename = 'eros3mill'      % str, name of high res model
        filenameLR = 'eros200700'   % str, name of low res model
        m = 6.687E15                % double, mass (kg)
        mu = G*m                    % double, gravitational parameter 
        alb = 0.8                   % double, albedo constant
        tEpoch = 2451170.5          % JD of epoch (days)
                                        % 1998 Dec 23 00:00:00 UTC
        oeHCI = [1.458260038106518*AU;
                 0.2228858603247133;
                 deg2rad(10.83015266864554);
                 deg2rad(304.4308844737856);
                 deg2rad(178.6132327246327);
                 deg2rad(208.1235381788443)];
                                    % 6x1 double, Heleocentric ecliptic OEs
                                        % [a e i Om w M_0]'
                                        % (km - rad rad rad rad)
        rotParam = [11.35;
                    17.22;
                    326.07;
                    1639.38864745];
                                    % 4x1 double, rotational parameters
                                        % [alpha delta prime rotRate]
                                        % (deg deg deg deg/day)
                                        % From IAU
    end
    
    methods
        function ast = ASTEROID()

        end
    end
    
end