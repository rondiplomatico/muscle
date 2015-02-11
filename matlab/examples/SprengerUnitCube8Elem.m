classdef SprengerUnitCube8Elem < muscle.AModelConfig
% Different tests for comparison between CMISS and KerMor.
%
% Variant 1: Fix the inner 5 nodes on each face of
    
    properties(SetAccess=private)
        Variant;
    end
    
    methods
        function this = SprengerUnitCube8Elem(varargin)
            varargin(end+1:end+2) = {'FL',2};
            this = this@muscle.AModelConfig(varargin{:});
            this.init;
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            m.T = 1;
            m.dt = .01;
            mu = m.DefaultMu;
            mu([1 2]) = [1 .1];
            m.DefaultMu = mu;
            
            %% Material configuration from CMISS/3Elem_sprenger.xml
            % c3M = 3.05907e-10 MPa
            mu(5) = 3.05907e-7; % [kPa]
            % c4M 
            mu(6) = 47.270456264135881; % [-]
            
            % c01 ^= c1M = 3.56463903963e-02 MPa
            mu(9) = 35.6463903963; % [kPa]
            % c01 ^= c2M = 3.859558659683e-3 MPa
            mu(10) = 3.859558659683; % [kPa]
            
            % sigma_max_calculation = 0.3 in MPa
            mu(13) = 300; % [kPa]
            % Exponential FL-fun
            mu(14) = .57; % [-]
            
            m.DefaultMu = mu;
        end
    end
    
    methods(Access=protected)
             
        function geo = getGeometry(this)
            s = load(fullfile(fileparts(which(mfilename)),'..','CMISS','Sprenger8Elem.mat'));
            geo = s.geo27;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            displ_dir(:,geo.Elements([1 3 5 7],geo.MasterFaces(1,:))) = true;
            displ_dir(:,geo.Elements([2 4 6 8],geo.MasterFaces(2,:))) = true;
        end
        
        function anull = seta0(~, anull)
            % Direction is xz
            anull(1,:,:) = 1;
        end
    end
    
    methods(Static)
        function test_SprengerUnitCube8Elem_ForceComparison
            % .1        .2          .3          .4
            % 14.9477   28.5314     41.0593     52.7679
            c = SprengerUnitCube8Elem;
            m = c.createModel;
            dfpos = m.getPositionDirichletBCFaceIdx([1 3 5 7],1);
            dfpos2 = m.getPositionDirichletBCFaceIdx([2 4 6 8],2);
            
            alphas = [.1 .2];% .3 .4];
            vis = [1 1];% 4 5];
            forces = zeros(2,length(alphas));
            pi = ProcessIndicator('Running tests for %d alpha values',length(alphas),false,length(alphas));
            mu = m.DefaultMu;
            for k = 1:length(alphas)
                mu(1) = vis(k);
                c.ActivationRampMax = alphas(k);
                [t,y] = m.simulate(mu);
                df = m.getResidualForces(t,y);
                forces(1,k) = sum(df(dfpos,end));
                forces(2,k) = sum(df(dfpos2,end));
                pi.step;
            end
            pi.stop;
            disp('Units are mN (milliNewton)');
            disp(forces);
        end
    end
    
end

