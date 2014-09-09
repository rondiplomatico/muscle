classdef SubMeshModelConfig < muscle.AModelConfig
    %SUBMESHMODELCONFIG A Model config providing the same settings as
    %another model config but for a submesh
    
    properties
        
    end
    
    properties(SetAccess=private)
        FullConfig;
        
        SubToFullNodeIdx;
        SubToFullElemIdx;
        SubToFullFaceIdx;
    end
    
    methods
        function this = SubMeshModelConfig(full, elems, varargin)
            geo = full.PosFE.Geometry;
            [subgeo, nodeidx, faceidx] = geo.getSubMesh(elems, varargin{:});
            this = this@muscle.AModelConfig(subgeo);
            
            % Store original config
            this.FullConfig = full;
            % Store index of elems/nodes in original geometry
            this.SubToFullNodeIdx = nodeidx;
            this.SubToFullElemIdx = elems;
            this.SubToFullFaceIdx = faceidx;
            
            % Copy other quantities
            this.NeumannCoordinateSystem = full.NeumannCoordinateSystem;
            this.a0CoordinateSystem = full.a0CoordinateSystem;
            this.FibreTypeWeights = full.FibreTypeWeights;
            this.Pool = full.Pool;
        end
        
        function configureModel(this, model)
            this.FullConfig.configureModel(model);
        end
        
        function prepareSimulation(this, mu, inputidx)
            this.FullConfig.prepareSimulation(mu, inputidx);
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            fullelemidx = this.SubToFullElemIdx(elemidx);
            fullfaceidx = this.SubToFullFaceIdx(faceidx);
            P = this.FullConfig.getBoundaryPressure(fullelemidx, fullfaceidx);
        end
        
        function u = getInputs(this)
            u = this.FullConfig.getInputs;
        end
        
        function x0 = getX0(this, x0)
            % todo
            x0 = this.FullConfig.getX0(x0);
        end
    end
    
    methods(Access=protected)
        function anull = seta0(this, ~)
            fullanull = this.FullConfig.geta0;
            anull = fullanull(:,:,this.SubToFullElemIdx);
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, ~, ~)
            [~, velo_dir, velo_dir_val] = this.FullConfig.getBC;
            velo_dir = velo_dir(:,this.SubToFullNodeIdx);
            velo_dir_val = velo_dir_val(:,this.SubToFullNodeIdx);
        end        
        
        function displ_dir = setPositionDirichletBC(this, ~)
            displ_dir = this.FullConfig.getBC;
            displ_dir = displ_dir(:,this.SubToFullNodeIdx);
        end
    end
    
end

