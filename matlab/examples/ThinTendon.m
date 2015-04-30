classdef ThinTendon < experiments.AExperimentModelConfig
    
    properties(SetAccess=private)
        ylen = 200;
        movetime = 80;
        configs;
    end
    
    methods
        function this = ThinTendon(varargin)
            this = this@experiments.AExperimentModelConfig(varargin{:});
            this.init;
            this.VelocityBCTimeFun = tools.ConstantUntil(this.movetime);
            this.configs = (-.9:.1:.2)*this.ylen/this.movetime;
            this.NumConfigurations = length(this.configs);
            this.NumOutputs = 1;
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            m.T = this.movetime+20;
            m.dt = m.T/200;
            
            mu = m.DefaultMu;
            % Small viscosity
            mu(3) = 0;
            m.DefaultMu = mu;
        end
        
        function tmr = getTendonMuscleRatio(this, points)
            % Pure tendon material!
            tmr = ones(1,size(points,2));            
        end
        
        function u = getInputs(this)
            % Returns the inputs `u(t)` of the model.
            %
            % if neumann boundary conditions are used, this input is
            % multiplied with the mu(3) parameter, which determines the
            % maximum force that is applied. u(t) determines its temporal
            % strength.
            %
            % this.Model can be used to get access to the model this
            % configuration is applied to.
            %
            % Return values:
            % u: The cell array of input functions to use within this
            % model.
            %
            % @type cell @default {}
            u = {this.getAlphaRamp(this.movetime,1)};
        end
        
        function P = getBoundaryPressure(~, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is megaPascal [MPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            if any(elemidx == 21:24) && faceidx == 4
                P = 1;
            end
        end
        
        function o = getOutputOfInterest(this, t, y)
            m = this.Model;
            % Get the residual dirichlet forces at end time
            df = m.getResidualForces(t(end),y(:,end));
            idx = m.getPositionDirichletBCFaceIdx(1:4,3);
            o(1) = sum(df(idx));
            %o(2) = sum(df(idx,passive_pos));   
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            geo = Belly.getBelly(6,this.ylen,'Radius',0,'InnerRadius',2,'Gamma',2);
            n = geo.Nodes;
            % Slightly deviate a node in the center
            n(1,137) = n(1,137)*1.1;
            geo = geometry.Cube27Node(n,geo.Elements);
            %geo = belly.scale(20);
            
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            % Fix ends in xz direction
            displ_dir(:,geo.Elements(1:4,geo.MasterFaces(3,:))) = true;
            displ_dir([1 3],geo.Elements(21:24,geo.MasterFaces(4,:))) = true;
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            % Fix ends in xz direction
            velo_dir(2,geo.Elements(21:24,geo.MasterFaces(4,:))) = true;
            velo_dir_val(velo_dir) = this.configs(this.CurrentConfigNr);
        end
    end
    
    methods(Static)
        function test_ThinTendon
            c = ThinTendon;
            m = c.createModel;
            c.CurrentConfigNr = 4;
            m.simulateAndPlot;
        end
        
        function runExperiment
            c = ThinTendon;
            m = c.createModel;
            e = tools.ExperimentRunner(m);
            d = e.runExperimentsCached;
            semilogy(d.o);
            title('Forces at fixed nodes at end time');
        end
    end
    
end

