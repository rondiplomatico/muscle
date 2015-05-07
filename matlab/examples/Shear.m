classdef Shear < muscle.AModelConfig
    % An example illustrating shear forces
    
    methods
        function this = Shear(varargin)
            this = this@muscle.AModelConfig(varargin{:});
            this.addOption('BC',1);
            this.init;
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            m.T = 50;
            m.dt = 1;
            m.ODESolver.RelTol = 1e-6;
            m.ODESolver.AbsTol = 1e-6;
            if this.Options.BC == 2
                m.DefaultInput = 1;
                m.DefaultMu(3) = 1;
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
            o = this.Options;
            if o.BC == 2
                if faceidx == 2 && ((o.GeoNr == 1 && elemidx == 1) || ...
                        (o.GeoNr == 2 && any(elemidx == 5:8)))
                   P = zeros(3,3);
                   P(3,1) = 1;
                end
            end
        end
        
        function u = getInputs(this)
            u = {};
            if this.Options.BC == 2
                u = {this.getAlphaRamp(this.Model.T/2,.1)};
            end
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            % Single cube with same config as reference element
            if this.Options.GeoNr == 1
                [pts, cubes] = geometry.Cube8Node.DemoGrid([0 1],[0 1],[0 1]);
            else
                [pts, cubes] = geometry.Cube8Node.DemoGrid(0:2,0:2,0:2);
            end
            geo = geometry.Cube8Node(pts, cubes);
            geo = geo.toCube27Node;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.PosFE.Geometry;
            o = this.Options;
            if o.GeoNr == 1
                displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
                if o.BC == 1
                    displ_dir(1:2,geo.Elements(1,geo.MasterFaces(2,:))) = true;
                end
            else
                displ_dir(:,geo.Elements(1:4,geo.MasterFaces(1,:))) = true;
                if o.BC == 1
                    displ_dir(1:2,geo.Elements(5:8,geo.MasterFaces(2,:))) = true;
                end
            end
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            % Determines the dirichlet velocities.
            %
            % The unit for the applied quantities is [mm/ms] = [m/s]
            geo = this.PosFE.Geometry;
            o = this.Options;
            if o.BC == 1
                if o.GeoNr == 1
                    velo_dir(3,geo.Elements(1,geo.MasterFaces(2,:))) = true;
                else
                    velo_dir(3,geo.Elements(5:8,geo.MasterFaces(2,:))) = true;
                end
                velo_dir_val(velo_dir) = .1;
            end
        end
        
        function anull = seta0(~, anull)
            anull(3,:,:) = 1;
        end
    end
    
    methods(Static)
        function test_Shear
            for g = 1:2
                for b = 1:2
                    f = Shear('GeoNr',g,'BC',b);
                    m = f.createModel;
                    m.plotGeometrySetup;
                    [t,y] = m.simulate;
                    m.plot(t,y,'Lambdas',true);
                end
            end
        end
    end
    
end

