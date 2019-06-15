classdef SPACECRAFT
    % Spacecraft object specifying orbit and camera properties
    %{
    -----------------------------------------------------------------------
    DESCRIPTION:
        Creates an object to store all of the ground truth data and
        simulation parameters for the ANS project. 
    -----------------------------------------------------------------------
    REFERENCES:
        - Text
    -----------------------------------------------------------------------
    NOTES:
        - Text
    -----------------------------------------------------------------------
    AUTHOR: Kaitlin Dennison
    -----------------------------------------------------------------------
    COPYTRIGHT: 2019 SLAB Group
    -----------------------------------------------------------------------
    TIMESTAMPS:
        - 09-Jun-2019: creation (KD)
        - 13-Jun-2019: style changes (KD)
        - 14-Jun-2019: style update and added properties (KD)
    -----------------------------------------------------------------------
    %}
    
    properties (SetAccess = private)
        name                % str, name of the sc 
        camera              % camera, short-range camera used by the sc
        tEpoch              % double, Julian date of epoch (days)
        oeACI               % 6x1 double, Asteroid Centered Inertial orbital 
                            % elements (navigation frame)
                                % [a e i Om w M_0]' 
                                % (km - rad rad rad rad)
        m                   % double, mass (kg)
        A                   % double, effective area (m^2)
        Cr                  % double, reflectivity coefficient
        FMR                 % double, fuel-mass ratio
    end
    methods
        function sc = SPACECRAFT(name,tEpoch,oeACI,camera)
            % Spacecraft object constructor
            %{
            ---------------------------------------------------------------
            INPUT:
                name:       str, name of test or simulation to add to any 
                            files created
                tEpoch:     double, Julian date of epoch (days)
                oeACI:      6x1 double, Asteroid Centered Inertial orbital 
                            elements (navigation frame)
                                [a e i Om w M_0]' 
                                (km - rad rad rad rad)
                camera:     camera, camera object associated with the 
                            spacecraft (See CAMERA.m)
            ---------------------------------------------------------------
            OUTPUT:
                sc:         SPACECRAFT object
            ---------------------------------------------------------------
            %}
            sc.name = name;
            sc.tEpoch = tEpoch;
            sc.oeACI = oeACI;
            sc.camera = camera;
        end
    end
end