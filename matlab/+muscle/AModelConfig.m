classdef AModelConfig < handle
    %AModelConfig
    
    properties
        Model;
    end
    
    properties(SetAccess=private)
        PosFE;
        
        PressFE;
    end
    
    properties(SetAccess=protected)
        FibreTypeWeights = [];
        Pool;
    end
    
    methods
        function this = AModelConfig(geo)
            if isa(geo,'geometry.Cube8Node')
                pos_geo = geo.toCube27Node;
                press_geo = geo;
            elseif isa(geo,'geometry.Cube20Node') || isa(geo,'geometry.Cube27Node')
                pos_geo = geo;
                press_geo = geo.toCube8Node;
            else
                error('Scenario not yet implemented for geometry class "%s"', class(geo));
            end
            if isa(pos_geo,'geometry.Cube27Node')
                this.PosFE = fem.HexahedronTriquadratic(pos_geo);
            else
                this.PosFE = fem.HexahedronSerendipity(pos_geo);
            end
            this.PressFE = fem.HexahedronTrilinear(press_geo);
        end
        
        function configureModel(this, model)
            % Overload this method to set model-specific quantities like
            % simulation time etc
            
            % do nothing by default
        end
        
        function P = getBoundaryPressure(~, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
        end
        
        function u = getInputFunction(this, m)
            u = @(t)1;
        end
        
        function x0 = getX0(this, x0)
            %% do nothing
        end
    end
    
    methods(Access=protected)
        function anull = seta0(~, anull)
            % do nothing!
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(~, velo_dir, velo_dir_val)
            % Determines the dirichlet velocities.
            %
            % The unit for the applied quantities is [mm/ms] = [m/s]
            %
            % In the default implementation there are no velocity
            % conditions.
        end
        
        function alpha = getAlphaRamp(~, ramptime, alphamax)
            if nargin < 3
                alphamax = 1;
            end
            alpha = @(t)alphamax * ((t<ramptime)*t/ramptime + (t>=ramptime));
        end
    end
    
    methods(Abstract, Access=protected)
        displ_dir = setPositionDirichletBC(this, displ_dir);
    end
    
    methods(Sealed)
        function [displ_dir, velo_dir, velo_dir_val] = getBC(this)
            N = this.PosFE.Geometry.NumNodes;
            displ_dir = false(3,N);
            displ_dir = this.setPositionDirichletBC(displ_dir);
            velo_dir = false(3,N);
            velo_dir_val = zeros(3,N);
            [velo_dir, velo_dir_val] = this.setVelocityDirichletBC(velo_dir, velo_dir_val);
            
            if any(any(displ_dir & velo_dir)) 
                error('Cannot impose displacement and velocity dirichlet conditions on same DoF');
            end
        end
        
        function anull = geta0(this)
            fe = this.PosFE;
            g = fe.Geometry;
            anull = zeros(3,fe.GaussPointsPerElem,g.NumElements);
            anull = this.seta0(anull);
            % Normalize anull vectors
            for m = 1:g.NumElements
                anull(:,:,m) = anull(:,:,m) ./ ([1;1;1]*Norm.L2(anull(:,:,m)));
            end
        end
    end
    
end

