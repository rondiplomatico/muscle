classdef Belly < muscle.AModelConfig
    
    properties(Access=private)
        NumParts;
    end
    
    methods
        function this = Belly
            np = 4;
            belly = Belly.getBelly(np,10,1,.5,2);
            this = this@muscle.AModelConfig(belly);
            this.NumParts = np;
            
            %% Muscle fibre weights
            types = [0 .2 .4 .6 .8 1];
            ftw = zeros(this.PosFE.GaussPointsPerElem,length(types),this.PosFE.Geometry.NumElements);
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
        end
        
        function configureModel(~, model)
            model.T = 150;
            model.dt = .1;
            f = model.System.f;
            f.alpha = @(t)0;
        end
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            for k = geo.NumElements-3:geo.NumElements
                displ_dir(:,geo.Elements(k,geo.MasterFaces(4,:))) = true;
            end
            for k = 1:4
                displ_dir(:,geo.Elements(k,geo.MasterFaces(3,:))) = true;
            end
        end
        
%         function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
%             %% Dirichlet conditions: Position (fix one side)
%             geo = this.PosFE.Geometry;
%             for k = [1 2 7 8]  
%                 velo_dir(1,geo.Elements(k,[1:3 9 10 13:15])) = true;
%             end
%             velo_dir_val(velo_dir) = -.1;
%         end
        
        function anull = seta0(~, anull)
            anull(2,:,:) = 1;
        end
    end
    
    methods(Static)
        function geo27 = getBelly(parts, length, radius, inner_radius, gamma)
            if nargin < 5
                gamma = 2;
                if nargin < 4
                    inner_radius = .5;
                    if nargin < 3
                        radius = 1;
                         if nargin < 2
                             length = 10;
                             if nargin < 1
                                 parts = 4;
                             end
                         end
                    end
                end
            end
            k = kernels.GaussKernel;
            k.Gamma = gamma;
            x = linspace(-length/2,length/2,parts*2+1);
            fx = k.evaluate(x,0)*radius;

            c = sin(pi/4);
            c8 = cos(pi/8);
            s8 = sin(pi/8);
            baseplane(:,:,1) = [0 .5 1  0 .5*c c8 0 s8 c
                                0  0 0 .5 .5*c s8 1 c8 c];
            R = [0 -1; 1 0];
            baseplane(:,:,2) = R*baseplane(:,[7 4 1 8 5 2 9 6 3],1);
            baseplane(:,:,3) = R*baseplane(:,[7 4 1 8 5 2 9 6 3],2);
            baseplane(:,:,4) = R*baseplane(:,[7 4 1 8 5 2 9 6 3],3);
            
            npp = 27*parts;
            nodes = zeros(3,npp*4);
            for p = 1:parts
                elempos = (p-1)*2 + (1:3);
                r = ones(9,1) * (inner_radius + fx(elempos)');
                z = ones(9,1) * x(elempos);
                for rotpart = 1:4
                    nodepos = ((p-1)*4 + (rotpart-1))*27 + (1:27);
                    nodes(:,nodepos) = [bsxfun(@times,repmat(baseplane(:,:,rotpart),1,3),r(:)'); z(:)'];
                end
            end

            % Filter nodes to unique ones
            [nodes, ~, elems] = unique(nodes','rows','stable');
            nodes = nodes';
            elems = reshape(elems,27,[])';
            geo27 = geometry.Cube27Node(nodes,elems);
            geo27.swapYZ;
        end
    end
    
end

