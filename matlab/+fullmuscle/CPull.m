classdef CPull < fullmuscle.AModelConfig
    
    properties(SetAccess=private)
        Version;
    end
    
    methods
        function this = CPull(version, xpts, ypts)
            % Creates a Debug simple muscle model configuration.
            %
            % Single cube with same config as reference element
            if nargin < 3
                ypts = [0 10];
                if nargin < 2
                    if nargin < 1
                        version = 1;
                    end
                    xpts = [0 10];
                    if version == 3
                        xpts = [0 10 20];
                    end
                end            
            end
            [pts, cubes] = geometry.Cube8Node.DemoGrid(xpts,ypts,[0 10]);
            geo = geometry.Cube8Node(pts, cubes);
            this = this@fullmuscle.AModelConfig(geo.toCube27Node);
            this.NeumannCoordinateSystem = 'global';
            this.Version = version;
        end
        
        function configureModel(this, m)
            configureModel@fullmuscle.AModelConfig(this, m);
            
            switch this.Version
                case [1 2]
                    m.T = 100;
                    m.dt = .5;
                case 3
                    m.T = 1000;
                    m.dt = 1;
            end
            
            m.DefaultMu = [1; 0; 1; 0];
            
            m.System.f.Pmax = 250;
            m.EnableTrajectoryCaching = true;
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            elem = 1;
            if this.Version == 3
                elem = 2;
            end
            if elemidx == elem && faceidx == 2
                P = 1;
            end
        end
        
        function u = getInputs(this)
            u{1,1} = this.getAlphaRamp(10,1);
            u{1,2} = this.getAlphaRamp(100,1);
            u{1,3} = this.getAlphaRamp(300,1);
            u{1,4} = this.getAlphaRamp(10,1,200);
            u{1,5} = this.getAlphaRamp(100,1,200);
            u{1,6} = this.getAlphaRamp(300,1,200);
            u{1,7} = this.getAlphaRamp(1000,1);
            null = @(t)ones(size(t));
            u(2,1:7) = repmat({null},1,7);
        end
        
    end
    
    methods(Access=protected)
        
        function ft = getFibreTypes(this)
            switch this.Version
                case 1
                   ft = 0;
                case 2
                   ft = [0 1];
                case 3
                   ft = [0 .2 .4 .6 .8 1];
            end
        end
        
        function sp = getSpindlePos(this)
            % Spindle position: first row element, second row gauss point
            % within element
            switch this.Version
                case 1
                   sp = [1; 1];
                case 2
                   sp = [1 1; 1 2];
                case 3
                   sp = [1 1 1 2 2 2; 1 10 20 4 10 25]; 
            end
        end
        
        function ftw = getFibreTypeWeights(this)
            % Get pre-initialized all zero weights
            ftw = getFibreTypeWeights@fullmuscle.AModelConfig(this);

            switch this.Version
                case 1
                   ftw(:,1,:) = 1;
                case 2
                   ftw(:,1,:) = .5;
                   ftw(:,2,:) = .5;
                case 3
                   fac = exp((1:6)/2);
                   fac = fac / sum(fac);
                   ftw(:,1,:) = fac(1);
                   ftw(:,2,:) = fac(2);
                   ftw(:,3,:) = fac(3);
                   ftw(:,4,:) = fac(4);
                   ftw(:,5,:) = fac(5);
                   ftw(:,6,:) = fac(6);
            end
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.PosFE.Geometry;
            % Always fix back side
            displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
        end
        
%         function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
%             % Determines the dirichlet velocities.
%             %
%             % The unit for the applied quantities is [mm/ms] = [m/s]
%             geo = this.PosFE.Geometry;
%             switch this.Version
%             
%             case {4,5,6}
%                 % Move the whole front
%                 velo_dir(1,geo.Elements(1,geo.MasterFaces(2,:))) = true;
%                 velo_dir_val(velo_dir) = .05;
%             end
%         end
        
        function anull = seta0(~, anull)
            anull(1,:,:) = 1;
        end
    end
    
end

