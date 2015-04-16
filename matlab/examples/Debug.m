classdef Debug < muscle.AModelConfig
    % A simple configuration for Debug purposes.
    % 
    % Uses a single undeformed cube with triquadratic position shape
    % functions and trilinear shape functions for pressure.
    %
    %
    % Version 1: No fibres, and one fixed side. Simple "does it simulate"
    % test
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
    %
    % Version 13: Unit tests. Pull at a unit cube with 1MPa, so should give
    % force of 1N
    % Version 14: Unit tests. Fix both sides and activate to max of 1MPa.
    % Should give a force of 1N on both sides
    
    methods
        function this = Debug(varargin)
            % Quick fix for old-style DebugConfig call passing in the
            % version number directly
            if nargin == 1 && isscalar(varargin{1})
                varargin = {'Version',varargin{1}};
            end
            % Creates a Debug simple muscle model configuration.
            this = this@muscle.AModelConfig(varargin{:});
            this.addOption('Version',3);
            this.init;
            
            switch this.Options.Version
                case {4,5,6}
                    this.VelocityBCTimeFun = tools.ConstantUntil(10);
                    this.Options.FL = 2;
                case {7,8,9}
                    this.VelocityBCTimeFun = tools.ConstantUntil(10);
                case 14
                    % Use classical FL function with fl(1)=1
                    this.Options.FL = 3;
            end
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            sys = m.System;
            f = sys.f;
            m.T = 100;
            m.dt = 1;
            m.ODESolver.RelTol = .00001;
            m.ODESolver.AbsTol = .002;
            switch this.Options.Version
            case 1
                m.DefaultMu(2) = -1;
            case 2
                m.DefaultMu(1:2) = [50; .5];
            case 3
                m.T = 400;
                m.dt = 1;
                types = [0 .2 .4 .6 .8 1];
                fe = this.PosFE;
                geo = fe.Geometry;
                ftw = zeros(fe.GaussPointsPerElem,length(types),geo.NumElements);
                % Test: Use only slow-twitch muscles
                ftw(:,1,:) = .4;
                ftw(:,2,:) = .05;
                ftw(:,3,:) = .05;
                ftw(:,4,:) = .1;
                ftw(:,5,:) = .2;
                ftw(:,6,:) = .2;
                this.FibreTypeWeights = ftw;
                p = models.motorunit.Pool;
                p.FibreTypes = types;
                this.Pool = p;
                m.ODESolver.RelTol = .01;
                m.ODESolver.AbsTol = .1;
                m.Plotter.DefaultArgs = {'Pool',true};
                m.DefaultMu(4) = 4;
            case {4,5,6}
                %% Material configuration from CMISS/3Elem_sprenger.xml                
                % c3M = 3.05907e-10 MPa
                m.DefaultMu(5) = 3.05907e-10; % b1 [MPa]
                % c4M = 47.270456264135881
                m.DefaultMu(6) = 47.270456264135881; % d1 [-]
                % c1M = 3.56463903963e-02 MPa
                m.DefaultMu(9) = 3.56463903963e1; % c10 [MPa]
                % c2M = 3.859558659683e-3 MPa
                m.DefaultMu(10) = 3.859558659683; % c01 [MPa]
                
                % sigma_max_calculation 
                m.DefaultMu(13) = .300; % [MPa]
                % lambda_ofl_calculation
                m.DefaultMu(14) = 1.4; % [-]
                m.T = 20;
                m.dt = .5;
            case {7,8,9}
                %% Using the Material configuration of Heidlauf
                m.DefaultMu(5) = 0.00355439810963035e-3; % b1 [MPa]
                m.DefaultMu(6) = 12.660539325481963; % d1 [-]
                m.DefaultMu(9) = 6.352e-13; % c10 [MPa]
                m.DefaultMu(10) = 3.627e-3; % c01 [MPa]
                
                m.T = 50;
                m.dt = .5;
            case 10
                f.alpha = @(t)0;
                m.ODESolver.RelTol = .001;
                m.ODESolver.AbsTol = .03;
                m.T = 1;
                m.dt = .01;
                m.DefaultMu(1:3) = [0;0;1];
                m.DefaultInput = 1;
            case 11
                m.ODESolver.RelTol = .001;
                m.ODESolver.AbsTol = .05;
                m.T = 10;
                m.dt = .05;
            case 12
                m.DefaultMu(1:2) = [.01; .04];
                m.System.UseCrossFibreStiffness = true;
            case 14
                m.DefaultMu(2) = 10;
                m.DefaultMu(13) = 1;
                m.T = 20;
                m.dt = .05;
            end
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is megaPascal [MPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            if this.Options.Version == 10
%                 if any(elemidx - (9:16) == 0) && faceidx == 2
                if elemidx == 1 && faceidx == 2
                   P = 2;
                end
            end
        end
        
        function u = getInputs(this)
            u = {};
            if this.Options.Version == 10
                u = {this.getAlphaRamp(this.Model.T/2,1)};
            end
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid([0 1],[0 1],[0 1]);
            %[pts, cubes] = geometry.Cube8Node.DemoGrid([0:.5:1],[0:.5:1],[0:.5:1]);
            if this.Options.Version == 11 || this.Options.Version == 12
                [pts, cubes] = geometry.Cube8Node.DemoGrid([-1 1],[-1 1],[-1 1]);
                pts(1,[6 8]) = 2;
                theta = .3;
                R = [cos(theta) -sin(theta) 0
                     sin(theta) cos(theta)  0
                     0          0           1];
%                 R = [cos(theta) -sin(theta) 0
%                        0 1 0
%                      sin(theta) 0 cos(theta)];

                pts = circshift(R,[1 1])*R*pts;
            end
            geo = geometry.Cube8Node(pts, cubes);
            geo = geo.toCube27Node;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.PosFE.Geometry;
            % Always fix back side
            switch this.Options.Version
                case 10
                    displ_dir(1,geo.Elements(1,geo.MasterFaces(1,:))) = true;
                case 14
                    displ_dir(:,geo.Elements(1,geo.MasterFaces(1:2,:))) = true;
                otherwise
                    displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
                    %displ_dir(:,geo.Elements(1:4,geo.MasterFaces(1,:))) = true;
            end
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            % Determines the dirichlet velocities.
            %
            % The unit for the applied quantities is [mm/ms] = [m/s]
            geo = this.PosFE.Geometry;
            switch this.Options.Version
            
            case {4,5,6}
                % Move the whole front
                velo_dir(1,geo.Elements(1,geo.MasterFaces(2,:))) = true;
                %velo_dir(1,geo.Elements(5:8,geo.MasterFaces(2,:))) = true;
                velo_dir_val(velo_dir) = .05;
            case {7,8,9}
                % Move the whole front
                velo_dir(1,geo.Elements(1,geo.MasterFaces(2,:))) = true;
                velo_dir_val(velo_dir) = .04;
            case 11
                velo_dir(:,geo.Elements(1,geo.MasterFaces(2,:))) = true;
                hlp = velo_dir;
                hlp([1 3],:) = 0;
                velo_dir_val(hlp) = .04;
            end
        end
        
        function anull = seta0(this, anull)
            switch this.Options.Version
            case {2,3,5,8,12,14}
                anull(1,:,:) = 1;
            case {4,7}
                % No fibres
            case {6,9}
                % Stretch perpendicular to fibre direction
                anull(2,:,:) = 1;
            case {10,11}
                anull(:,:,:) = 0;
            end
        end
    end
    
    methods(Static)
        function test_DebugConfigV4_9(version)
            % Results:
            % 4: Force at T=20, x-position (center node) 1.49844: -0.0476511N
            % 5: Force at T=20, x-position (center node) 1.4986: -0.170556N
            % 6: Force at T=20, x-position (center node) 1.49844: -0.0476511N
            % 7: Force at T=50, x-position (center node) 1.39931: 0.00210816N
            % with  4^3 gaussp:  
            %       5^3 gaussp: 
            % 8: Force at T=50, x-position (center node) 1.39056: 0.000453817N
            % with  4^3 gaussp: 
            %       5^3 gaussp: 
            % 9: Force at T=50, x-position (center node) 1.39931: 0.00210815N
            % with 4^3 gaussp: 
            %      5^3 gaussp: 
            if nargin < 1
                version = 4;
            end
            m = muscle.Model(Debug(version));
            [t,y] = m.simulate;
            pdf = m.getResidualForces(t,y);
            
            % 1..27 are residual forces from fixed position nodes
            force_on_fixed_side = sum(pdf(1:27,:),1);
            % 28..36 (=9) are residual forces from moved nodes
            force_on_moved_side = sum(pdf(28:end,:),1);
            figure;
            plot(t,force_on_fixed_side,'r',t,force_on_moved_side,'b');
            xlabel('time'); ylabel('force [N]');

            % Node 15 is center -> dof indices (15-1)*3 + [1 2 3]
            yall = m.System.includeDirichletValues(t,y);
            xpos = yall((15-1)*3+1,end);
            fprintf('Force at T=%g, x-position (center node) %g: %gN\n',t(end),xpos,force_on_moved_side(end));
            
            m.plot(t,y,'DF',pdf,'Velo',true,'Pause',.02);
        end
        
        function test_Config13
            m = muscle.Model(Debug(13));
            [t,y] = m.simulate;
            pdf = m.getResidualForces(t,y);
            
            % 1..27 are residual forces from fixed position nodes
            force_on_fixed_side = sum(pdf(1:27,:),1);
            % 28..36 (=9) are residual forces from moved nodes
            force_on_moved_side = sum(pdf(28:end,:),1);
            figure;
            plot(t,force_on_fixed_side,'r',t,force_on_moved_side,'b');
            xlabel('time'); ylabel('force [N]');

            % Node 15 is center -> dof indices (15-1)*3 + [1 2 3]
            yall = m.System.includeDirichletValues(t,y);
            xpos = yall((15-1)*3+1,end);
            fprintf('Force at T=%g, x-position (center node) %g: %gN\n',t(end),xpos,force_on_moved_side(end));
            
            m.plot(t,y,'DF',pdf,'Velo',true,'Pause',.02);
        end
        
        function test_Config14
            m = muscle.Model(Debug(14));
            [t,y] = m.simulate;
            pdf = m.getResidualForces(t,y);
            
            % 1..27 are residual forces from fixed position nodes
            force_on_fixed_side = sum(pdf(1:27,:),1);
            % 28..36 (=9) are residual forces from moved nodes
            force_on_moved_side = sum(pdf(28:end,:),1);
            figure;
            plot(t,force_on_fixed_side,'r',t,force_on_moved_side,'b');
            xlabel('time'); ylabel('force [N]');

            % Node 15 is center -> dof indices (15-1)*3 + [1 2 3]
            yall = m.System.includeDirichletValues(t,y);
            xpos = yall((15-1)*3+1,end);
            fprintf('Force at T=%g, x-position (center node) %g: %gN\n',t(end),xpos,force_on_moved_side(end));
            
            m.plot(t,y,'DF',pdf,'Velo',true,'Pause',.02);
        end
        
        function test_DebugConfigs(versions)
            if nargin < 1
                versions = [1 2 3 10 11];
            end
            for k = versions
                fprintf('Testing DebugConfig Version %d\n',k);
                m = muscle.Model(Debug(k));
                [t,y] = m.simulate;
                df = m.getResidualForces(t,y);
                m.plot(t,y,'DF',df,'Velo',true);
            end
        end
    end
    
end

