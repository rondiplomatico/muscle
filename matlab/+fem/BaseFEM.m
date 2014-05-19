classdef BaseFEM < handle
    %FEMBASE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Geometry;
        
        % Values of basis functions on gauss points
        Ngp;
        
        % The mass matrix
        M;
                
        elem_detjac;
        face_detjac;
        Ngpface;
        dN_facenormals;
        NormalsOnFaceGP;
        
        % transformed basis function gradients, stored in eldofs x gp*3
        % x element matrix (accessible in 20x3-chunks for each gauss point and element)
        transgrad;
    end
        
    methods
        function this = BaseFEM(geometry)
            this.Geometry = geometry;
            this.init;
        end
        
        function init(this)
            g = this.Geometry;
            el = g.Elements;
            np = g.NumNodes;
            
            %% Precomputations
            % Iterate each cube
            gp = g.gaussp;
            % Number of gauss points
            ngp = g.GaussPointsPerElem;
            Mass = zeros(np,np);
            % Number of elements
            nel = size(el,1);
            % nodes per Element
            eldofs = size(el,2);
            
            % Basis function evaluation on Gauss points
            Ngpval = zeros(eldofs,ngp,nel);
            % Jacobian of element deformation at gauss points
            eljac = zeros(nel,ngp);
            % transformed basis function gradients, stored in eldofs x gp*3
            % x element matrix (accessible in 20x3-chunks for each gauss point and element)
            dtnall = zeros(eldofs,ngp*3,nel);
            % Iterate all volumes
            for m = 1:nel
                bval = zeros(eldofs,eldofs);
                elem = el(m,:);
                % Iterate all gauss points
                for gi = 1:ngp
                    % gauss point \xi
                    xi = gp(:,gi);
                    
                    %% N on gauss points
                    Ngpval(:,gi,m) = this.N(xi);
                    
                    %% Jacobian related
                    dNxi = this.gradN(xi);
                    % Jacobian of isogeometric mapping
                    Jac = g.Nodes(:,elem)*dNxi;
                    
                    %% Mass matrix related
                    % Jacobian at gauss point: get coordinates of nodes (8
                    % for trilinear, 20 for triquadratic) and multiply with
                    % gradient
                    eljac(m,gi) = det(Jac);
                    % Evaluate basis functions at xi
                    Nxi = this.N(xi);
                    % Add up [basis function values at xi] times [volume
                    % ratio at xi] times [gauss weight for xi]
                    bval = bval + g.gaussw(gi)*(Nxi*Nxi')*eljac(m,gi);
                    
                    %% Transformed Basis function gradients at xi
                    dtnall(:,3*(gi-1)+1:3*gi,m) = dNxi / Jac;
                end
                % Build up upper right part of symmetric mass matrix
                for j=1:eldofs
                    Mass(elem(j),elem(j:end)) = Mass(elem(j),elem(j:end)) + bval(j,j:end);
                end
            end
            this.Ngp = Ngpval;
            this.transgrad = dtnall;
            this.M = sparse(Mass + Mass' - diag(diag(Mass)));
            this.elem_detjac = eljac;
            
            %% Boundary/Face precomputations
            nf = g.NumFaces;
            ngp = g.GaussPointsPerElemFace;
            npf = g.NodesPerFace;
            Ngpval = zeros(npf,ngp,nf);
            facejac = zeros(nf,ngp);
            dNN = zeros(npf,ngp,nf);
            transNormals = zeros(3,ngp,nf);
            for fn = 1:nf
                elemidx = g.Faces(1,fn);
                faceidx = g.Faces(2,fn);
                masterfacenodeidx = g.MasterFaces(faceidx,:);
                % Get N values on all gauss points of the face, returning a
                % DofsPerElem x ngp vector.
                Nall = this.N(g.facegaussp(:,:,faceidx));
                % From that, we only need the values on the face nodes
                Ngpval(:,:,fn) = Nall(masterfacenodeidx,:);
                
                dNxi = this.gradN(g.facegaussp(:,:,faceidx));    
                for gi = 1:ngp
                    dNxipos = [0 9 18]+gi;
%                     facenodeidx = g.Elements(elemidx,masterfacenodeidx);
                    
                    % Get full transformation jacobian
                    Jac = g.Nodes(:,g.Elements(elemidx,:))*dNxi(:,dNxipos);

                    % Precompute transformed normals
                    transNormals(:,gi,fn) = Jac*g.FaceNormals(:,faceidx);
%                     transNormals(:,gi,fn) = g.FaceNormals(:,faceidx);
                    
                    % Take only the restriction of the jacobian to the
                    % face-dimensions for transformation theorem
                    Jac = Jac(g.FaceDims(:,faceidx),g.FaceDims(:,faceidx));
%                     Jac_check = g.Nodes(:,g.Elements(elemidx,:)) * this.gradN(xi);
%                     Jac_check = Jac_check(g.FaceDims(:,faceidx),g.FaceDims(:,faceidx));
%                     det(Jac) - det(Jac_check)
                    facejac(fn,gi) = det(Jac);
                    
                    % Precompute transformation of face normals, only the
                    % gradient N * n part (will be left-multiplied with u
                    % at the corresponding locations to get the deformed
                    % normal, can be used for plotting)
                    dNN(:,gi,fn) = dNxi(masterfacenodeidx,dNxipos) * g.FaceNormals(:,faceidx);
                end
                transNormals(:,:,fn) = Norm.normalizeL2(transNormals(:,:,fn));
            end
            this.dN_facenormals = dNN;
            this.face_detjac = facejac;
            this.Ngpface = Ngpval;
            this.NormalsOnFaceGP = transNormals;
        end
                
        function plot(this, pm)
            if nargin < 2
                pm = PlotManager;%(false,1,2);
                pm.LeaveOpen = true;
            end
            
            h = pm.nextPlot('mass','Mass matrix','node','node');
            mesh(h,full(this.M));
            
            pm.nextPlot('mass_pattern','Mass matrix pattern','dof','dof');
            spy(this.M);
            
            if nargin < 2
                pm.done;
            end
        end
    end
    
    methods(Abstract)
        % Evaluates the elementary basis functions on the geometry master
        % element
        %
        % Parameters:
        % x: A 3xn matrix with points `x_i` in columns @type matrix<double>
        %
        % Return values:
        % nx: The values `N_i(x_j)` of each elementary basis function `N_i` per
        % row on all given points `x_i` @type matrix<double>
        nx = N(this, x);
        
        % Evaluates the gradients of the elementary basis functions on the
        % master element
        %
        % Parameters:
        % x: A 3xn matrix with points `x_i` in columns @type matrix<double>
        %
        % Return values:
        % nx: The values `[\frac{\partial N}{x_1} N_i(x_j), \frac{\partial
        % N}{x_3} N_i(x_j), \frac{\partial N}{x_3} N_i(x_j)]` of each
        % elementary basis function `N_i` per row on all given points
        % `x_j`. If `x` is only a 3x1 column vector, this corresponds to
        % `\nabla N_i(x)` @type matrix<double>
        dnx = gradN(this, x);
    end
    
    methods(Static,Access=protected)
        function res = test_BasisFun(subclass)
            h = 1e-8;
            res = true;
            for k=1:10
                a = rand(3,1);
                A = repmat(a,1,3);
                diff = (subclass.N(A + h*eye(3)) - subclass.N(A))/h - subclass.gradN(a);
                res = res & all(diff(:) < 10*h);
            end
        end
    end
        
    methods(Static)
        function res = test_JacobiansDefaultGeo
            % Tests if the deformation jacobians using linear and quadratic
            % elements is the same for a default geometry
            ranges = {-1:1, -2:1, -2:2};
            res = true;
            for k = 1:length(ranges)
                [pts, cubes] = geometry.Cube8Node.DemoGrid(ranges{k},ranges{k});
                g = geometry.Cube8Node(pts, cubes);
                tl = fem.HexahedronTrilinear(g);
                tq = fem.HexahedronTriquadratic(g.toCube27Node);
                res = res && norm(tl.elem_detjac-tq.elem_detjac,'inf') < 1e-14;
            end
        end
    end
    
end

