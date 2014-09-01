classdef LongMORexample < muscle.AModelConfig
% A long geometry with 20% deviation from default cubic positions and
% complex fibre structure

    properties(Constant)
       OutputDir = fullfile(fileparts(which(mfilename)),'LongMOR_output') 
    end
    %
    %%
    properties  % TODO: set access to what?
        
        % element divisions in x-, y-, z-direction
        Partsx;
        Partsy;
        Partsz;

    end
    %%
    methods
        function this = LongMORexample(devi)
            if nargin < 1
                devi = .2;
            end
            px = -10:2.5:10;
            py = -40:2.5:40;
            pz = 0:2.5:10;
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube20Node.DemoGrid(px, py, pz, devi);
            geo = geometry.Cube20Node(pts, cubes);
            this = this@muscle.AModelConfig(geo);
            this.Partsx = px;
            this.Partsy = py;
            this.Partsz = pz;
        end
        
        function configureModel(this, model)
            %
            model.T = 100;
            model.dt = 1;
            model.EnableTrajectoryCaching = true;
            
            model.Data.useFileTrajectoryData;
            model.ComputeTrajectoryFxiData = true;
            
            % default model parameters 
            model.DefaultMu = [1; 50];      % mu = [viscosity; activation duration]
            
            % specify model parameters (mu = [viscosity; activation duration])
            sys = model.System;
            %viscosity
            sys.Params(1).Range = [0.5 10];
            sys.Params(1).Desired = 10;
            sys.Params(1).Spacing = 'log';
            % activation duration
            sys.Params(2).Name = 'alpha-ramp';
            sys.Params(2).Range = [20 100];
            sys.Params(2).Desired = 10;
            sys.Params(2).Spacing = 'log';
            
            f = sys.f;
            % Material set (see main comment)
            f.c10 = 6.352e-10; % [kPa]
            f.c01 = 3.627; % [kPa]
            f.b1 = 2.756e-5; % [kPa]
            f.d1 = 43.373; % [-]

            f.Pmax = 73; % [kPa]
            f.lambdafopt = 1.2; % [-]            
        end
        
        function prepareSimulation(this, mu, inputidx)
            % 
            sys = this.Model.System;
            sys.f.alpha = this.getAlphaRamp(mu(2),1);    % (in ..ms, up to maxvalue.., starting at ..ms)
        end
                
    end
    %
    %%
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            % Dirichlet conditions: Position (fix both sides)
            geo = this.PosFE.Geometry;
            
            NrEltsx = size(this.Partsx,2) -1;
            NrEltsy = size(this.Partsy,2) -1;
            NrEltsz = size(this.Partsz,2) -1;            
            
            for ix = 1:NrEltsx
                for iz = 1:NrEltsz
                    % left side
                    displ_dir(:,geo.Elements((ix-1)*NrEltsy*NrEltsz+iz,geo.MasterFaces(3,:))) = true;
                    % right side
                    displ_dir(:,geo.Elements(ix*NrEltsy*NrEltsz-(iz-1),geo.MasterFaces(4,:))) = true;
                end
            end

        end
        
        function anull = seta0(this, anull)
            fe = this.PosFE;
            if fe.GaussPointsPerElem ~= 27
                warning('a0 designed for 27 gauss points!');
            end
            
            NrEltsx = size(this.Partsx,2) -1;
            NrEltsy = size(this.Partsy,2) -1;
            NrEltsz = size(this.Partsz,2) -1;
            
            x = linspace(0,1,NrEltsy);
            basea0 = [-sin(x*pi); cos(x*pi)];
            for ix = 1:NrEltsx
                for iy = 1:NrEltsy
                    for iz = 1:NrEltsz
                        anull(2,:,(ix-1)*NrEltsy*NrEltsz+(iy-1)*NrEltsz+iz) = basea0(2,iy);
                        anull(3,:,(ix-1)*NrEltsy*NrEltsz+(iy-1)*NrEltsz+iz) = basea0(1,iy);
                    end
                end
            end
            
        end
        
    end
    
end