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
        function this = LongMORexample(varargin)
            % specify as a subclass of AModelConfig
            this = this@muscle.AModelConfig(varargin{:});
            % initialise -- calls getGeometry()
            this.init;
        end
        
        function configureModel(this, model)
            % vary tolerances
            model.ODESolver.RelTol = 1e-3;
            model.ODESolver.AbsTol = 1e-3;
            
            % simulation time and step size
            model.T = 100;
            model.dt = 1;
            
            % for offline phases
            model.EnableTrajectoryCaching = true;
            
            model.Data.useFileTrajectoryData;
            model.ComputeTrajectoryFxiData = true;
            
            % change default model parameters mu
            model.DefaultMu(1) = 1;             % viscosity
            model.DefaultMu(2) = 50;            % activation duration
            model.DefaultMu(3) = 0;             % NeumannBC (max force) = default
            model.DefaultMu(9) = 6.352e-10;     % Mooney-Rivlin c10
            model.DefaultMu(10) = 3.627;        % Mooney-Rivlin c01
            model.DefaultMu(13) = 73;           % muscle fibre maximal force
            
            % specify model parameters (mu = [viscosity; activation duration])
            sys = model.System;
            % viscosity
            sys.Params(1).Range = [0.5 10];
            sys.Params(1).Desired = 10;
            sys.Params(1).Spacing = 'log';
            % activation duration
            sys.Params(2).Range = [20 100];
            sys.Params(2).Desired = 10;
            sys.Params(2).Spacing = 'log';
            
            % simulation: in mu(2) ms up to maxvalue 1, starting at t=0
            this.ActivationRampMax = 1;     % default
            this.ActivationRampOffset = 0;  % default
        end
                
        function [oy,oz] = getOutputOfInterest(this, model, t, uvw)
            % Writes the data of interest into oy and oz
            %
            geo = model.Config.PosFE.Geometry;
            %[df,nf] = model.getResidualForces(t,uvw);
            uvw = model.System.includeDirichletValues(t,uvw);
            %
            facenode_idx = [];
            % get displacement (x,y,z) of middle node (element 445, face 5, node 8)
            facenode_idx = [facenode_idx; model.getFaceDofsGlobal(445,5)];
            % get displacement (x,y,z) of node at 1/3rd (element 413, face 5, node 8)
            facenode_idx = [facenode_idx; model.getFaceDofsGlobal(413,5)];
            %
            o = uvw(facenode_idx,:);            % gives matrix, where each rowvector shows x-,y- or z-positon of one node over time
            % middle node, track z-direction
            oz = o(24,:);
            % 1/3rd node, track y-direction
            oy = o(47,:);
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
        
        function geo = getGeometry(this)
            % geometry for model configuration
            px = -10:2.5:10;
            py = -40:2.5:40;
            pz = 0:2.5:10;
            devi = .2;
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube20Node.DemoGrid(px, py, pz, devi);
            geo = geometry.Cube20Node(pts, cubes);
            
            this.Partsx = px;
            this.Partsy = py;
            this.Partsz = pz;
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