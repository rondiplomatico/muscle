classdef BaseFEM < handle
    %FEMBASE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Dependent)
        GaussPointRule;
        GaussPointsPerElem;
        GaussPointsPerElemFace;
    end
    
    properties(SetAccess=private)
        Geometry;
        
        % Values of basis functions on gauss points
        Ngp;
        
        % The mass matrix
        M;
        
        % The damping matrix (for linear damping)
        D;
        
        % The element-dof to global node assembly matrix
        Sigma;
                
        elem_detjac;
        face_detjac;
        Ngpface;
        dN_facenormals;
        
        % transformed basis function gradients, stored in eldofs x gp*3
        % x element matrix (accessible in 20x3-chunks for each gauss point and element)
        transgrad;
        
        %% Gauss integration related properties
        GaussPoints;
        GaussWeights;
        FaceGaussPoints;
        FaceGaussWeights;
        NormalsOnFaceGP;
        FaceAreas;
    end
    
    properties(Access=private)
        fGaussPointRule = 3;
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
            
            % Get the right gauss points
            this.initGaussPoints;
            
            % Build element-dof to global dof assembly matrix
            flat_el = el';
            flat_el = flat_el(:)';
            i = []; j = [];
            for k=1:np
                idx = find(flat_el == k);
                i = [i ones(size(idx))*k];%#ok
                j = [j idx];%#ok
            end
            this.Sigma = sparse(i,j,ones(size(i)),np,g.NumElements*g.DofsPerElement);
            
            %% Precomputations
            % Number of gauss points
            ngp = this.GaussPointsPerElem;
            Mass = zeros(np,np);
            Damp = Mass;
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
                mass_gp = zeros(eldofs,eldofs);
                damp_gp = zeros(eldofs,eldofs);
                elem = el(m,:);
                % Iterate all gauss points
                for gi = 1:ngp
                    % gauss point \xi
                    xi = this.GaussPoints(:,gi);
                    
                    %% N on gauss points
                    Ngpval(:,gi,m) = this.N(xi);
                    
                    %% Jacobian related
                    dNxi = this.gradN(xi);
                    % Jacobian of isogeometric mapping
                    Jac = g.Nodes(:,elem)*dNxi;
                    
                    %% Mass matrix related
                    % Jacobian at gauss point
                    eljac(m,gi) = abs(det(Jac));
                    
                    % Evaluate basis functions at xi
                    Nxi = this.N(xi);
                    
                    % Add up [basis function values at xi] times [volume
                    % ratio at xi] times [gauss weight for xi]
                    mass_gp = mass_gp + this.GaussWeights(gi)*(Nxi*Nxi')*eljac(m,gi);
                    
                    % Add up damping values
                    damp_gp = damp_gp + this.GaussWeights(gi)*(dNxi*dNxi')*eljac(m,gi);
                    
                    %% Transformed Basis function gradients at xi
                    dtnall(:,3*(gi-1)+1:3*gi,m) = dNxi / Jac;
                end
                % Build up upper right part of symmetric mass matrix
                for j=1:eldofs
                    Mass(elem(j),elem(j:end)) = Mass(elem(j),elem(j:end)) + mass_gp(j,j:end);
                    Damp(elem(j),elem(j:end)) = Damp(elem(j),elem(j:end)) + damp_gp(j,j:end);
                end
            end
            this.Ngp = Ngpval;
            this.transgrad = dtnall;
            this.M = sparse(Mass + Mass' - diag(diag(Mass)));
            this.D = sparse(Damp + Damp' - diag(diag(Damp)));
            this.elem_detjac = eljac;
            
            %% Boundary/Face precomputations
            nf = g.NumFaces;
            ngp = this.GaussPointsPerElemFace;
            npf = g.NodesPerFace;
            Ngpval = zeros(npf,ngp,nf);
            facejac = zeros(nf,ngp);
            dNN = zeros(npf,ngp,nf);
            transNormals = zeros(3,ngp,nf);
            % tangent dimensions (1st+2nd row) and normal dimension index
            % (3rd row)
            tangidx = [2 2 1 1 1 1;
                       3 3 3 3 2 2];
            signum = [-1 1 1 -1 -1 1];
            for fn = 1:nf
                elemidx = g.Faces(1,fn);
                faceidx = g.Faces(2,fn);
                masterfacenodeidx = g.MasterFaces(faceidx,:);
                % Get N values on all gauss points of the face, returning a
                % DofsPerElem x ngp vector.
                Nall = this.N(this.FaceGaussPoints(:,:,faceidx));
                % From that, we only need the values on the face nodes
                Ngpval(:,:,fn) = Nall(masterfacenodeidx,:);
                
                dNxi = this.gradN(this.FaceGaussPoints(:,:,faceidx));    
                for gi = 1:ngp
                    dNxipos = [0 ngp 2*ngp]+gi;
                    
                    % Get full transformation jacobian
                    Jac = g.Nodes(:,g.Elements(elemidx,:))*dNxi(:,dNxipos);

                    % Get rotational part of jacobian
                    [q, ~] = qr(Jac);
                    
                    % Precompute transformed normals
                    % sign(Jac(tangidx(3,fn),tangidx(3,fn)))...
                    transNormals(:,gi,fn) = signum(faceidx)...
                        *cross(Jac(:,tangidx(1,faceidx)),Jac(:,tangidx(2,faceidx)));
                    
                    % Take out rotational part (stretch only)
                    Jac = q\Jac;
                    % Take only the restriction of the jacobian to the
                    % face-dimensions for transformation theorem
                    Jac = Jac(g.FaceDims(:,faceidx),g.FaceDims(:,faceidx));
                    
                    facejac(fn,gi) = abs(det(Jac));
                    
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
            this.FaceAreas = facejac*this.FaceGaussWeights;
            this.Ngpface = Ngpval;
            this.NormalsOnFaceGP = transNormals;
        end
        
        function a = getFaceArea(this, elemidx, faceidx)
            g = this.Geometry;
            a = 0;
            for k = 1:length(elemidx)
                a = a + this.FaceAreas(g.Faces(1,:) == elemidx(k) & g.Faces(2,:) == faceidx(k));
            end
        end
        
        function v = getElementVolume(this, elemidx)
            % Returns the volume of the element with the specified index in
            % [mm³]
            v = this.elem_detjac(elemidx,:)*this.GaussWeights;
        end
        
        function v = getTotalVolume(this)
            % Returns the total geometry volume in [mm³]
            v = sum(this.elem_detjac*this.GaussWeights);
        end
        
        function gp = getGlobalGaussPoints(this, elemidx)
            % Returns the positions of all gauss points of the specified
            % element in global (=reference) coordinates
            g = this.Geometry;
            gp = g.Nodes(:,g.Elements(elemidx,:)) * this.Ngp(:,:,elemidx);
        end
                
        function plot(this, pm)
            if nargin < 2
                pm = PlotManager(false,2,2);
                pm.LeaveOpen = true;
            end
            
            h = pm.nextPlot('mass','Mass matrix','node','node');
            mesh(h,full(this.M));
            
            pm.nextPlot('mass_pattern','Mass matrix pattern','dof','dof');
            spy(this.M);
            
            h = pm.nextPlot('damp','Damping matrix','node','node');
            mesh(h,full(this.D));
            
            pm.nextPlot('damp_pattern','Damping matrix pattern','dof','dof');
            spy(this.D);
            
            if nargin < 2
                pm.done;
            end
        end
        
        function value = get.GaussPointRule(this)
            value = this.fGaussPointRule;
        end
        
        function set.GaussPointRule(this, value)
            if ~isequal(this.fGaussPointRule, value)
                
                this.fGaussPointRule = value;
                this.init;
            end
        end
        
        function nc = get.GaussPointsPerElem(this)
            nc = size(this.GaussPoints,2);
        end
        
        function nc = get.GaussPointsPerElemFace(this)
            nc = size(this.FaceGaussPoints,2);
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
    
    methods(Access=private)
        function initGaussPoints(this)
            %% Init Gauss points
            switch this.fGaussPointRule
                case 3 % 3-point-rule            
                    g = [-sqrt(3/5) 0 sqrt(3/5)];
                    w = [5/9 8/9 5/9];
                case 4 % 4-point-rule
                    g = [-sqrt(3/7 + 2/7*sqrt(6/5)) -sqrt(3/7 - 2/7*sqrt(6/5)) sqrt(3/7 - 2/7*sqrt(6/5)) sqrt(3/7 + 2/7*sqrt(6/5))];
                    w = [(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36];
                case 5 % 5-point rule
                    a = 2*sqrt(10/7);
                    g = [-sqrt(5+a)/3 -sqrt(5-a)/3 0 sqrt(5-a)/3 sqrt(5+a)/3];
                    a = 13*sqrt(70);
                    w = [(322-a)/900 (322+a)/900 128/225 (322+a)/900 (322-a)/900];
                case 6
                    % values from http://pomax.github.io/bezierinfo/legendre-gauss.html
                    v = [0.3607615730481386	0.6612093864662645
                        0.3607615730481386	-0.6612093864662645
                        0.4679139345726910	-0.2386191860831969
                        0.4679139345726910	0.2386191860831969
                        0.1713244923791704	-0.9324695142031521
                        0.1713244923791704	0.9324695142031521];
                    g = v(:,2)';
                    w = v(:,1)';
                case 7
                    v = [0.4179591836734694	0.0000000000000000
                        0.3818300505051189	0.4058451513773972
                        0.3818300505051189	-0.4058451513773972
                        0.2797053914892766	-0.7415311855993945
                        0.2797053914892766	0.7415311855993945
                        0.1294849661688697	-0.9491079123427585
                        0.1294849661688697	0.9491079123427585];
                    g = v(:,2)';
                    w = v(:,1)';
                case 8
                    v = [0.3626837833783620	-0.1834346424956498
                    0.3626837833783620	0.1834346424956498
                    0.3137066458778873	-0.5255324099163290
                    0.3137066458778873	0.5255324099163290
                    0.2223810344533745	-0.7966664774136267
                    0.2223810344533745	0.7966664774136267
                    0.1012285362903763	-0.9602898564975363
                    0.1012285362903763	0.9602898564975363];
                    g = v(:,2)';
                    w = v(:,1)';
                case 9
                    v = [0.3302393550012598	0.0000000000000000
                    0.1806481606948574	-0.8360311073266358
                    0.1806481606948574	0.8360311073266358
                    0.0812743883615744	-0.9681602395076261
                    0.0812743883615744	0.9681602395076261
                    0.3123470770400029	-0.3242534234038089
                    0.3123470770400029	0.3242534234038089
                    0.2606106964029354	-0.6133714327005904
                    0.2606106964029354	0.6133714327005904];
                    g = v(:,2)';
                    w = v(:,1)';
                case 10
                    v = [0.2955242247147529	-0.1488743389816312
                    0.2955242247147529	0.1488743389816312
                    0.2692667193099963	-0.4333953941292472
                    0.2692667193099963	0.4333953941292472
                    0.2190863625159820	-0.6794095682990244
                    0.2190863625159820	0.6794095682990244
                    0.1494513491505806	-0.8650633666889845
                    0.1494513491505806	0.8650633666889845
                    0.0666713443086881	-0.9739065285171717
                    0.0666713443086881	0.9739065285171717];
                    g = v(:,2)';
                    w = v(:,1)';
                case 15
                    v = [0.2025782419255613	0.0000000000000000
                    0.1984314853271116	-0.2011940939974345
                    0.1984314853271116	0.2011940939974345
                    0.1861610000155622	-0.3941513470775634
                    0.1861610000155622	0.3941513470775634
                    0.1662692058169939	-0.5709721726085388
                    0.1662692058169939	0.5709721726085388
                    0.1395706779261543	-0.7244177313601701
                    0.1395706779261543	0.7244177313601701
                    0.1071592204671719	-0.8482065834104272
                    0.1071592204671719	0.8482065834104272
                    0.0703660474881081	-0.9372733924007060
                    0.0703660474881081	0.9372733924007060
                    0.0307532419961173	-0.9879925180204854
                    0.0307532419961173	0.9879925180204854];
                    g = v(:,2)';
                    w = v(:,1)';
                case 20
                    v = [0.1527533871307258	-0.0765265211334973
                    0.1527533871307258	0.0765265211334973
                    0.1491729864726037	-0.2277858511416451
                    0.1491729864726037	0.2277858511416451
                    0.1420961093183820	-0.3737060887154195
                    0.1420961093183820	0.3737060887154195
                    0.1316886384491766	-0.5108670019508271
                    0.1316886384491766	0.5108670019508271
                    0.1181945319615184	-0.6360536807265150
                    0.1181945319615184	0.6360536807265150
                    0.1019301198172404	-0.7463319064601508
                    0.1019301198172404	0.7463319064601508
                    0.0832767415767048	-0.8391169718222188
                    0.0832767415767048	0.8391169718222188
                    0.0626720483341091	-0.9122344282513259
                    0.0626720483341091	0.9122344282513259
                    0.0406014298003869	-0.9639719272779138
                    0.0406014298003869	0.9639719272779138
                    0.0176140071391521	-0.9931285991850949
                    0.0176140071391521	0.9931285991850949];
                    g = v(:,2)';
                    w = v(:,1)';
                otherwise
                    error('Quadrature rule for %d points per dimension not implemented\nYou can obtain more coefficients from e.g. http://pomax.github.io/bezierinfo/legendre-gauss.html',this.fGaussPointRule);
            end
            
            %% Transfer to 3D
            [WX,WY,WZ] = meshgrid(w);
            [GX,GY,GZ] = meshgrid(g);
            W = WX.*WY.*WZ;
            this.GaussPoints = [GX(:) GY(:) GZ(:)]';
            this.GaussWeights = W(:);
            
            % Faces Gauss points
            [WX,WY] = meshgrid(w);
            [GX,GY] = meshgrid(g);
            W = WX.*WY;
            fgp = zeros(3,length(g)^2,6);
            faceremainingdim = [-1 1 -1 1 -1 1];
            g = this.Geometry;
            for faceidx = 1:6
                fgp(g.FaceDims(:,faceidx),:,faceidx) = [GX(:) GY(:)]';
                fgp(~g.FaceDims(:,faceidx),:,faceidx) = faceremainingdim(faceidx);
            end
            this.FaceGaussPoints = fgp;
            this.FaceGaussWeights = W(:);
        end
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

