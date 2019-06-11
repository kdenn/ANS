classdef ANSGT
    % Ground truth associated with the entire ANS simulation
    % Kaitlin Dennison - SLAB - Spring 2019
    
    % Constructor Options:
        % asteroid: 
    
    properties (SetAccess = protected)
        % Simulation Parameters
        asteroid            % asteroid structure
        nSpacecraft         % number of spacecraft
        spacecraft          % cell array of spacecraft objects [nSpacecraftx1]
        tEpoch              % Julian date of epoch (days)
        tInterval           % time interval between images (s)
        nImages             % number of images per spacecraft
        % Ground Truth Data
        tJD                 % Julian date of each image (days) [nImagesx1]
        stateGT             % ground truth of state estimate [(nSpacecraftxN+M)xnImages]
        images              % images corresponding to the stateGT
        features            % features corresponding to each image
        % Feature D&T Parameters
        
    end
    methods
        function ansGT = ANSGT(varargin)
            % ANSGT constructor
            %% Defaults
            AU = 149597870.7; % AU conversion [km/AU]
            G = 6.67259*10^-20; % Gravitational Constant
            % Asteroid
                defaultAst.name = 'Eros';
                defaultAst.filename = 'eros3mill'; % High Res Model
                defaultAst.filenameLR = 'eros200700'; % Low Res Model 
                defaultAst.m = 6.687E15; % mass (kg)
                defaultAst.mu = G*defaultAst.m;
                defaultAst.alb = 0.8;  % albedo
                % Orbital Elements (JPL HORIZONS)
                defaultAst.JD_epoch = 2451170.5;   % JD epoch (days)
                    % 1998 Dec 23 00:00:00 UTC
                defaultAst.oeSCI = [1.458260038106518*AU;
                                    0.2228858603247133;
                                    deg2rad(10.83015266864554);
                                    deg2rad(304.4308844737856);
                                    deg2rad(178.6132327246327);
                                    deg2rad(0.5596952381222570)/(24*60*60);
                                    deg2rad(208.1235381788443)];
                    % [a e i Om w n M_0]' 
                    % (km - rad rad rad rad/s rad)
                % Rotation Parameters (IAU Report) epoch: J2000
                defaultAst.rotParam = [11.35;
                                       17.22;
                                       326.07;
                                       1639.38864745];
                	% [ascention, declination, prime meridian, rotation rate]' 
                    % (deg deg deg deg)
            % Spacecraft
                tEpoch = datenum(2019,9,4,0,0,0)+1721058.5;  % JD epoch (days)
                oeACI = [35 0.01 deg2rad(80) deg2rad(220) 0 1.5]';
                defaultSc = SPACECRAFT(tEpoch,oeACI,defaultAst.mu,'defaultANS');
                defaultNSc = 2;            
            % Simulation
                defaultTEp = tEpoch;
                defaultTInt = 5*60;
                defaultNIm = 500;
                
            %% Input Parser
            p = inputParser;
            addParameter(p,'asteroid',defaultAst)
            addParameter(p,'nSpacecraft',defaultNSc,@isnumeric)
            addParameter(p,'spacecraft',{})
            addParameter(p,'camera','default',@ischar)
            addParameter(p,'tEpoch',defaultTEp,@isnumeric)
            addParameter(p,'tInterval',defaultTInt,@isnumeric)
            addParameter(p,'nImages',defaultNIm,@isnumeric)
            parse(p,varargin{:})
            
            %% Check Input
            switch p.Results.camera
                case 'PtGrey17'
                    % Point Grey 17mm
                    defaultSc.camera = CAMERA(17,[1200 1920],[5.86 5.86],'PtGrey17');
                case 'NEAR'
                    % NEAR Shoemaker
                    defaultSc.camera = CAMERA(168,[244 537],[27 16],'NEAR');
                case 'OSIRIS'
                    % OSIRIS REx's NavCam
                    defaultSc.camera = CAMERA(7.6,[1944 2592],[2.2 2.2],'NavCam');
                case 'NanoCam'
                    % GOMspace NanoCam
                    defaultSc.camera = CAMERA(8,[2048 1536],[3.2 3.2],'NanoCam');
                otherwise
                    % Point Grey
                    defaultSc.camera = CAMERA(12,[1200 1920],[5.86 5.86],'PtGrey12');
            end
            if isempty(p.Results.spacecraft)
                sc = cell(p.Results.nSpacecraft,1);
                for n = 1:p.Results.nSpacecraft
                    sc{n} = defaultSc;
                end
                ansGT.spacecraft = sc;
                ansGT.nSpacecraft = p.Results.nSpacecraft;
                for n = 1:ansGT.nSpacecraft
                    ansGT.spacecraft{n}.name = ['orbiter',num2str(n)];
                end
            else
                ansGT.nSpacecraft = size(p.Results.spacecraft,1);
                ansGT.spacecraft = p.Results.spacecraft;
            end
            
            %% Parse Remaining Input
            ansGT.asteroid = p.Results.asteroid;
            ansGT.tEpoch = p.Results.tEpoch;
            ansGT.tInterval = p.Results.tInterval;
            ansGT.nImages = p.Results.nImages;
            tIntD = ansGT.tInterval/86400;
            ansGT.tJD = (ansGT.tEpoch:tIntD:(ansGT.tEpoch+tIntD*(ansGT.nImages-1)))';
        end
        
        function setProp(ansGT,varargin)
            if mod(nargin,2) ~= 0
                error('An even number of inputs must be given')
            end
            for n = 1:2:nargin
                switch varargin{n}
                    case 'asteroid'
                        ansGT.asteroid = varargin{n+1};
                    otherwise
                        disp(['input: ',varargin{n}, ' not recognized'])
                end
            end
            
        end
        
        function ansGT = compGT(ansGT)
            % Compute all of the ground truth meta data (not the images)
        end
        
        function ansGT = genIms(ansGT)
            % Generate all of the images
            % compGT must be called first!
        end
        
        function ansGT = detectFeat(ansGT)
            % Detect features in all of the images 
            % genIms must be called first!
        end
    end
end