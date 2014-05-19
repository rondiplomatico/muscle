classdef SprengerUnitCube8Elem < muscle.AModelConfig
% Different tests for comparison between CMISS and KerMor.
%
% Variant 1: Fix the inner 5 nodes on each face of
    
    properties(SetAccess=private)
        Variant;
    end
    
    methods
        function this = SprengerUnitCube8Elem(variant)
            if nargin < 1
                variant = 1;
            end
            s = load(fullfile(fileparts(which(mfilename)),'..','CMISS','Sprenger8Elem.mat'));
            this = this@muscle.AModelConfig(s.geo20,s.geo8);
            this.Variant = variant;
        end
        
        function configureModel(~, model)
            model.T = 1;
            model.dt = .01;
            f = model.System.f;
            f.Viscosity = 1;
            
            %% Material configuration from CMISS/3Elem_sprenger.xml
            % malpha_calculation
            f.alpha = 1; % [-]
            
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
        end
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            if this.Variant == 1
                % 1:4 one side
                displ_dir(:,geo.Elements(1,[11 16 18])) = true;
                displ_dir(:,geo.Elements(3,[9 13 16])) = true;
                displ_dir(:,geo.Elements(5,[4 6 11])) = true;
                displ_dir(:,geo.Elements(7,[1 4 9])) = true;
%                 displ_dir(1,geo.Elements(1,[1 4 6 9 13])) = true;
%                 displ_dir(1,geo.Elements(3,[1 4 6 11 18])) = true;
%                 displ_dir(1,geo.Elements(5,[1 9 13 16 18])) = true;
%                 displ_dir(1,geo.Elements(7,[6 11 13 16 18])) = true;

                % 29:32 other side
                displ_dir(:,geo.Elements(2,[12 17 20])) = true;
                displ_dir(:,geo.Elements(4,[10 15 17])) = true;
                displ_dir(:,geo.Elements(6,[5 8 12])) = true;
                displ_dir(:,geo.Elements(8,[3 5 10])) = true;                
%                 displ_dir(1,geo.Elements(2,[3 5 8 10 15])) = true;
%                 displ_dir(1,geo.Elements(4,[3 5 8 12 20])) = true;
%                 displ_dir(1,geo.Elements(6,[3 10 15 17 20])) = true;
%                 displ_dir(1,geo.Elements(8,[8 12 15 17 20])) = true;
            else
                displ_dir(:,geo.Elements(1,[1 4 6 9 11 13 16 18])) = true;
                displ_dir(:,geo.Elements(3,[1 4 6 9 11 13 16 18])) = true;
                displ_dir(:,geo.Elements(5,[1 4 6 9 11 13 16 18])) = true;
                displ_dir(:,geo.Elements(7,[1 4 6 9 11 13 16 18])) = true;
            end
            
        end
        
        function anull = seta0(~, anull)
            % Direction is xz
            anull(1,:,:) = 1;
        end
    end
    
    methods(Static)
        function test_ForceComparisonVariant1
            % .1        .2          .3          .4
            % 14.9477   28.5314     41.0593     52.7679
            m = muscle.Model(SprengerUnitCube8Elem(1));
            f = m.System.f;
            mu = 1;
            g = m.Config.PosFE.Geometry;
            % Get node indices from one side
            idx = unique([g.Elements(1,[11 16 18]) g.Elements(3,[9 13 16])...
                g.Elements(5,[4 6 11]) g.Elements(7,[1 4 9])]);
            % Compute DoF positions of nodes (x,y,z)
            dofidx = [(idx-1)*3+1 (idx-1)*3+2 (idx-1)*3+3];
            % Get relative position from overall dirichlet index
            [~, relpos] = intersect(m.System.bc_dir_displ_idx, dofidx);
            [~, otherrelpos] = setdiff(m.System.bc_dir_displ_idx, dofidx);
            
            alphas = [.1 .2];% .3 .4];
            vis = [1 1];% 4 5];
            forces = zeros(2,length(alphas));
            pi = ProcessIndicator('Running tests for %d alpha values',length(alphas),false,length(alphas));
            for k = 1:length(alphas)
                f.alpha = alphas(k);
                f.Viscosity = vis(k);
                [t,y] = m.simulate(mu);
                [pdf] = m.getResidualForces(t,y);
                forces(1,k) = sum(abs(pdf(relpos,end)));
                forces(2,k) = sum(abs(pdf(otherrelpos,end)));
                pi.step;
            end
            pi.stop;
            disp('Units are mN (milliNewton)');
            disp(forces);
        end
        
        function test_ForceComparisonVariant2
            % .1        .2          .3          .4
            % 14.9477   28.5314     41.0593     52.7679
            m = muscle.Model(SprengerUnitCube8Elem(2));
            f = m.System.f;
            mu = 1;
            % Select central node on loose side
            idx = 15;
            % Compute DoF positions of nodes (x,y,z)
            dofidx = [(idx-1)*3+1 (idx-1)*3+2 (idx-1)*3+3];
            
            
            alphas = [.1 .2];% .3 .4];
            vis = [1 1];% 4 5];
            pos = zeros(3,length(alphas));
            pi = ProcessIndicator('Running tests for %d alpha values',length(alphas),false,length(alphas));
            for k = 1:length(alphas)
                f.alpha = alphas(k);
                f.Viscosity = vis(k);
                [t,y] = m.simulate(mu);
                y = m.System.includeDirichletValues(y);
                pos(:,k) = y(dofidx,end);
                pi.step;
            end
            pi.stop;
            disp('Units are mN (milliNewton)');
            disp(forces);
        end
    end
    
end

