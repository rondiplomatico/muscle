classdef DebugConfig < fullmuscle.AModelConfig
    % A simple configuration for Debug purposes.
    % 
    % Uses a single undeformed cube with triquadratic position shape
    % functions and trilinear shape functions for pressure.
    %
    %
    % Version 1: No fibres, and one fixed side
    % Version 2: Fibres, but no FibreTypes
    % Version 3: Fibres with fibre type distribution and APExp
    % evaluations
    %
    % Versions 4-9 are related to one-element tests, where the element is
    % stretched in x-direction over time and the forces at
    % both the fixed and moved side are extracted.
    % Tests 4-6 use material parameters from M. Sprenger and stretch to
    % x=1.5[mm], and 7-9 use the built-in Material from T. Heidlauf and
    % stretch to 1.4[mm].
    % Version 4/7: No fibres
    % Version 5/8: Fibres in x-direction
    % Version 6/9: Fibres in y-direction
    %
    % Version 10: One side is fixed and a certain pressure is applied to
    % the other side's face
    % Version 11: A slightly rotated cube with one side fixed and velocity
    % boundary conditions
    % Version 12: A cube with fibres and additional cross-fibre direction
    % markert pressure. (only activation used)
    
    properties(SetAccess=private)
        Version;
    end
    
    methods
        function this = DebugConfig(version)
            % Creates a Debug simple muscle model configuration.
            %
            if nargin < 1
                version = 1;
            end
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid([0 1],[0 1],[0 1]);
            geo = geometry.Cube8Node(pts, cubes);
            this = this@fullmuscle.AModelConfig(geo.toCube27Node);
            
            this.Version = version;
        end
        
        function configureModel(this, m)
            configureModel@fullmuscle.AModelConfig(this, m);
            
            switch this.Version
            case 1
                m.T = 100;
                m.dt = .1;
                
                m.ODESolver.RelTol = .001;
                m.ODESolver.AbsTol = .01;
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
            if this.Version == 10
%                 if any(elemidx - (9:16) == 0) && faceidx == 2
                if elemidx == 1 && faceidx == 2
                   P = 2;
                end
            end
        end
    end
    
    methods(Access=protected)
        
        function ft = getFibreTypes(~)
%             ft = [0 .2 .4 .6 .8 1];
%             ft = [0 .4 1];
            ft = 0;
        end
        
        function ftw = getFibreTypeWeights(this)
            % Get pre-initialized all zero weights
            ftw = getFibreTypeWeights@fullmuscle.AModelConfig(this);

            ftw(:,1,:) = 1;
%             ftw(:,1,:) = .3;
%             ftw(:,2,:) = .4;
%             ftw(:,3,:) = .3;
%                 ftw(:,4,:) = .1;
%                 ftw(:,5,:) = .2;
%                 ftw(:,6,:) = .2;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.PosFE.Geometry;
            % Always fix back side
%             if this.Version == 10
%                 displ_dir(1,geo.Elements(1,geo.MasterFaces(1,:))) = true;
%             else
                displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
%                 displ_dir(:,geo.Elements(1,geo.MasterFaces(2,:))) = true;
%                 displ_dir(:,geo.Elements(1,geo.MasterFaces(3,:))) = true;
%             end
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
            switch this.Version
            case {1}
                anull(1,:,:) = 1;
%             case {4,7}
%                 % No fibres
%             case {5,8}
%                 % Stretch along fibre direction
%                 anull(1,:,:) = 1;
%             case {6,9}
%                 % Stretch perpendicular to fibre direction
%                 anull(2,:,:) = 1;
%             case {10,11}
%                 anull(:,:,:) = 0;
            end
        end
    end
    
    methods(Static)
        function test_DebugConfig(version)
            if nargin < 1
                version = 1;
            end
            m = muscle.Model(muscle.DebugConfig(version));
            [t,y] = m.simulate;
            df = m.getResidualForces(t,y);
            
            m.plot('DF',df);
        end
    end
    
end

