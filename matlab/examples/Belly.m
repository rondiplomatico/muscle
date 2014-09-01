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
        function geo27 = getBelly(parts, length, radius, inner_radius, gamma, minz)
            if nargin < 6 || minz == 0
                minz = [];
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
            end
            x = linspace(0,length,parts*2+1);
            if isa(radius,'function_handle')
                fx = radius(x);
                if size(fx,1) == 1
                    fx = [fx; fx];
                end
            else
                if isscalar(gamma)
                    gamma = [gamma gamma];
                end
                if isscalar(radius)
                    radius = [radius radius];
                end
                k = kernels.GaussKernel;
                k.Gamma = gamma(1);
                fx(1,:) = k.evaluate(x,length/2)'*radius(1);
                k.Gamma = gamma(2);
                fx(2,:) = k.evaluate(x,length/2)'*radius(2);
            end
            if isscalar(inner_radius)
                inner_radius = [1 1]*inner_radius;
            end

            c = sin(pi/4);
            c8 = cos(pi/8);
            s8 = sin(pi/8);
            baseplane(:,:,1) = [0 .5 1  0 .5*c c8 0 s8 c
                                0  0 0 .5 .5*c s8 1 c8 c];
            % Mirror first segment
            baseplane(:,:,2) = [-1 0; 0 1]*baseplane(:,[3 2 1 6 5 4 9 8 7],1);
            % Rotate others by 180°
            baseplane(:,:,3) = -baseplane(:,9:-1:1,1);
            baseplane(:,:,4) = -baseplane(:,9:-1:1,2);
            
            npp = 27*parts;
            nodes = zeros(3,npp*4);
            for p = 1:parts
                elempos = (p-1)*2 + (1:3);
                rx = ones(9,1) * (inner_radius(1) + fx(1,elempos));
                ry = ones(9,1) * (inner_radius(2) + fx(2,elempos));
                z = ones(9,1) * x(elempos);
                for rotpart = 1:4
                    nodepos = ((p-1)*4 + (rotpart-1))*27 + (1:27);
                    nodes(1,nodepos) = bsxfun(@times,repmat(baseplane(1,:,rotpart),1,3),rx(:)');
                    nodes(2,nodepos) = bsxfun(@times,repmat(baseplane(2,:,rotpart),1,3),ry(:)');
                    nodes(3,nodepos) = z(:);
                end
            end

            % Filter nodes to unique ones
            [nodes, ~, elems] = unique(nodes','rows','stable');
            nodes = nodes';
            elems = reshape(elems,27,[])';
            
            % This little hack allows to restrict the minimum z node
            % position to a certain value, so as if the muscle belly is
            % "lying on a table"
            %
            % This only works on the "downside" faces, so if the threshold
            % is larger than the intermediate node position of the side
            % faces, the muscle geometry will appear "dented".
            if ~isempty(minz)
                hlp_g27 = geometry.Cube27Node;
                centerline_facenodeidx = hlp_g27.MasterFaces(3,:);
                for elem = [3:4:4*parts 4:4:4*parts]
                    centerline_nodes = elems(elem,centerline_facenodeidx);
                    toScale = nodes(2,centerline_nodes) < minz;
                    factor = minz./nodes(2,centerline_nodes(toScale));
                    nodes(2,centerline_nodes(toScale)) = minz;
                    % Also scale inner node positions
                    inner_nodes = elems(elem,[4:6 13:15 22:24]);
                    nodes(2,inner_nodes(toScale)) = nodes(2,inner_nodes(toScale)).*factor;
                end
            end
            
            geo27 = geometry.Cube27Node(nodes,elems);
            geo27.swapYZ;
        end
        
        function res = test_BellyGeometrieGeneration
            g = Belly.getBelly;
            g = Belly.getBelly(4,35,1,.2,7);
            g = Belly.getBelly(4,35,1,[.2 .6],5);
            g = Belly.getBelly(4,35,1,.2,[10 20]);
            g = Belly.getBelly(4,35,1,[.2 .6],[10 20]);
            g = Belly.getBelly(4,35,[4 2],.5,[10 20]);
            g = Belly.getBelly(4,35,@(x)[sqrt(abs(x)); 1./(x-34).^2],.2);
            g = Belly.getBelly(4,35,@(x)sqrt(abs(x)));
            g.plot;
            g = Belly.getBelly(4,35,[4 2],.5,[10 20],-1.5);
            g.plot;
            res = 1;
        end
    end
    
end

