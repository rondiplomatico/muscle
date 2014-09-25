classdef CPull < fullmuscle.AModelConfig
    
    methods
        function this = CPull
            % Creates a Debug simple muscle model configuration.
            %
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid([0 1],[0 1],[0 1]);
            geo = geometry.Cube8Node(pts, cubes);
            this = this@fullmuscle.AModelConfig(geo.toCube27Node);
            this.NeumannCoordinateSystem = 'global';
        end
        
        function configureModel(this, m)
            configureModel@fullmuscle.AModelConfig(this, m);
            m.T = 100;
            m.dt = .1;
                        
            m.DefaultMu = [1; 0; 1; 0];
            
            m.System.f.Pmax = 250;
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            if elemidx == 1 && faceidx == 2
                P = 1;
            end
        end
        
        function u = getInputs(this)
            u{1} = this.getAlphaRamp(10,1);
            u{2} = this.getAlphaRamp(100,1);
            u{3} = this.getAlphaRamp(300,1);
        end
        
    end
    
    methods(Access=protected)
        
        function ft = getFibreTypes(~)
%             ft = [0 .2 .4 .6 .8 1];
%             ft = [0 .1];
            ft = 0;
        end
        
        function sp = getSpindlePos(~)
            % Spindle position: first row element, second row gauss point
            % within element
%             sp = [1 1; 1 2];
            sp = [1; 1];
        end
        
        function ftw = getFibreTypeWeights(this)
            % Get pre-initialized all zero weights
            ftw = getFibreTypeWeights@fullmuscle.AModelConfig(this);

            ftw(:,1,:) = 1;
%             ftw(:,2,:) = .5;
%             ftw(:,2,:) = .4;
%             ftw(:,3,:) = .3;
%                 ftw(:,4,:) = .1;
%                 ftw(:,5,:) = .2;
%                 ftw(:,6,:) = .2;
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
        
        function anull = seta0(this, anull)
                anull(1,:,:) = 1;
        end
    end
    
end

