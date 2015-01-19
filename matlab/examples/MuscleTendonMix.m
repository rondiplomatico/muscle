classdef MuscleTendonMix < muscle.AModelConfig
    % Muscle - Tendon mixed geometries
    %
    % Contains several examples:
    % Variant 1: Simple 
    
    properties
        Variant;
    end
    
    methods
        function this = MuscleTendonMix(variantnr)
            if nargin < 1
                variantnr = 1;
            end
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid(0:1,-1:3,0:1);
            geo = geometry.Cube8Node(pts, cubes);
            this = this@muscle.AModelConfig(geo.toCube20Node);
            this.Variant = variantnr;
        end
        
        function configureModel(this, m)
            m.T = 4;
            m.dt = .01;
            m.DefaultInput = 1;
            f = m.System.f;
            f.alpha = @(t)0;
            
            os = m.ODESolver;
            os.RelTol = .0001;
            os.AbsTol = .05;
            
            switch this.Variant
                case {1 2}
                    m.DefaultMu = [1; 0; 40; 0; 2.756e-5; 43.373; 7.99; 16.6];
            end
            
            % Ramp up the external pressure
            m.System.Inputs{1} = this.getAlphaRamp(1,1);
        end
        
        function tmr = getTendonMuscleRatio(this)
            % Returns the [0,1] ratio between tendon and muscle at all
            % gauss points of all elements
            %
            % This method simply returns an all-zero ratio, meaning muscle only. 
            tmr = getTendonMuscleRatio@muscle.AModelConfig(this);
            switch this.Variant
                case 1
                    tmr(:,1) = 0;%.4;
                    tmr(:,2) = .33;%.4;
                    tmr(:,3) = .66;%.4;
                    tmr(:,4) = 1;%.4;
                case 2 
                    tmr(19:27,:) = 1;%.4;
            end
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            switch this.Variant
                % Pull on front face
                case {1 2}
                    if elemidx == 1 && faceidx == 3
                        P = -1;
                    end
                %case 2
                %    if elemidx == 1 && faceidx == 1
                %        P = 1;
                %    end
            end
        end
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            
            switch this.Variant
                case {1 2}
                    % Fix all on left and only the y,z directions of the back right
                    % nodes
                    displ_dir(2,geo.Elements(4,geo.MasterFaces(4,:))) = true;
            end
            %displ_dir(:,geo.Elements(4,geo.MasterFaces(4,:))) = true;
            
            %displ_dir(1,geo.Elements(3,[9 18 27])) = false;
            %displ_dir(1,geo.Elements(2,[8 16 20])) = false;
        end
        
        function anull = seta0(~, anull)
           % Direction is y
           anull(2,:,:) = 1;
        end
    end
    
end

