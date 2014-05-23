classdef DebugConfig < muscle.AModelConfig
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
    
    properties(SetAccess=private)
        Version;
    end
    
    methods
        function this = DebugConfig(version)
            % Creates a Debug simple muscle model configuration.
            %
            if nargin < 1
                version = 3;
            end
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid([0 1],[0 1],[0 1]);
            geo = geometry.Cube8Node(pts, cubes);
            this = this@muscle.AModelConfig(geo.toCube27Node,geo);
            
            this.Version = version;
        end
        
        function configureModel(this, model)
            sys = model.System;
            f = sys.f;
            model.T = 1;
            model.dt = .05;
            model.ODESolver.RelTol = .00001;
            model.ODESolver.AbsTol = .002;
            switch this.Version
            case 1
                f.alpha = 0;
            case 3
                model.T = 400;
                model.dt = 1;
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
                p = motorunit.Pool;
                p.FibreTypes = types;
                this.Pool = p;
                model.ODESolver.RelTol = .01;
                model.ODESolver.AbsTol = .1;
            case {4,5,6}
                %% Material configuration from CMISS/3Elem_sprenger.xml
                % c1M = 3.56463903963e-02 MPa
                f.c10 = 35.6463903963; % [kPa]
                % c2M = 3.859558659683e-3 MPa
                f.c01 = 3.859558659683; % [kPa]
                % c3M = 3.05907e-10 MPa
                f.b1 = 3.05907e-7; % [kPa]
                % c4M 
                f.d1 = 47.270456264135881; % [-]
                % sigma_max_calculation = 0.3 in MPa
                f.Pmax = 300; % [kPa]
                % lambda_ofl_calculation
                f.lambdafopt = 1.4; % [-]
                f.ForceLengthFun = @(ratio)(ratio<=1).*exp(-((1-ratio)/.57).^4) + (ratio>1).*exp(-((ratio-1)/.14).^3);
                % The derivative of the force-length function as function handle
                f.ForceLengthFunDeriv = @(ratio)(ratio<=1).*((1/.57)*(((1-ratio)/.57).^3).*exp(-((1-ratio)/.57).^4)) ...
                    - (ratio > 1) .* ((1/.14) .* (((ratio-1)/.14).^2) .* exp(-((ratio-1)/.14).^3));
                model.T = 12;
                model.dt = .2;
                f.alpha = 0;
                sys.ApplyVelocityBCUntil = 10;
            case {7,8,9}
                %% Using the "normal" Material configuration of Heidlauf
                model.T = 15;
                model.dt = .1;
                f.alpha = 0;
                sys.ApplyVelocityBCUntil = 10;
            end
        end
        
    end
    
    methods(Access=protected)
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.PosFE.Geometry;
            % Always fix back side
            displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
            switch this.Version
            
            case {2,3}
                %displ_dir(:,geo.Elements(1,[6:8 11 12 18:20])) = true;
            
                % Only fix the right back corner nodes in yz-direction
    %             displ_dir([2 3],geo.Elements(1,[9 18 27])) = true;
%                 displ_dir(:,geo.Elements(1,[9 18 27])) = true;
            end
            
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            % Determines the dirichlet velocities.
            %
            % The unit for the applied quantities is [mm/ms] = [m/s]
            geo = this.PosFE.Geometry;
            switch this.Version
            
            case {4,5,6}
                % Move the whole front
                velo_dir(1,geo.Elements(1,geo.MasterFaces(2,:))) = true;
                velo_dir_val(velo_dir) = .05;
            case {7,8,9}
                % Move the whole front
                velo_dir(1,geo.Elements(1,geo.MasterFaces(2,:))) = true;
                velo_dir_val(velo_dir) = .04;
            end
        end
        
        function anull = seta0(this, anull)
            switch this.Version
            case {2,3}
                anull(1,:,:) = 1;
            case {4,7}
                % No fibres
            case {5,8}
                % Stretch along fibre direction
                anull(1,:,:) = 1;
            case {6,9}
                % Stretch perpendicular to fibre direction
                anull(2,:,:) = 1;
            end
        end
    end
    
    methods(Static)
        function test_DebugConfigV4_9(version)
            % Results:
            % 4: Force at T=12, x-position (center node) 1.5: -89.4587mN
            % 5: Force at T=12, x-position (center node) 1.5: -222.511mN
            % 6: Force at T=12, x-position (center node) 1.5: -89.458mN
            %
            % 7: Force at T=15, x-position (center node) 1.40088: -5.11973mN
            % with  4^3 gaussp:  -5.12226mN
            %       5^3 gaussp: -5.12224mN
            % 8: Force at T=15, x-position (center node) 1.4: -98.0233mN
            % with  4^3 gaussp:  -101.762mN
            %       5^3 gaussp: -102.297mN
            % 9: Force at T=15, x-position (center node) 1.4: -5.11482mN
            % with 4^3 gaussp:  -5.12227mN, 5^3 gaussp: -5.11725mN
            if nargin < 1
                version = 4;
            end
            m = muscle.Model(muscle.DebugConfig(version));
            mu = 1;
            [t,y] = m.simulate(mu);
            pdf = m.getResidualForces(t,y);
            
            % 1..27 are residual forces from fixed position nodes
            force_on_fixed_side = sum(pdf(1:27,:),1);
            % 28..36 (=9) are residual forces from moved nodes
            force_on_moved_side = sum(pdf(28:end,:),1);
            figure;
            plot(t,force_on_fixed_side,'r',t,force_on_moved_side,'b');
            xlabel('time'); ylabel('force [mN]');

            % Node 15 is center -> dof indices (15-1)*3 + [1 2 3]
            yall = m.System.includeDirichletValues(t,y);
            xpos = yall((15-1)*3+1,end);
            fprintf('Force at T=%g, x-position (center node) %g: %gmN\n',t(end),xpos,force_on_moved_side(end));
            
            %m.plot(t,y,'PDF',pdf,'Velo',true);
%             m.plot(t,y,'PDF',pdf);
        end
    end
    
end

