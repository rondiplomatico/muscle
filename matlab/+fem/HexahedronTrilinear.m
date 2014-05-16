classdef HexahedronTrilinear < fem.BaseFEM
    % HexahedronTrilinear: Base class for linear ansatz functions on
    % hexahedral geometry (8-Point elements)
    
    methods
        function this = HexahedronTrilinear(geo)
            if nargin < 1
                geo = geometry.Cube8Node;
            end
            this = this@fem.BaseFEM(geo);
            %this.EdgeIndices = 1:8;
        end
    
        function Nx = N(~, x)
            Nx = [(1-x(1,:)).*(1-x(2,:)).*(1-x(3,:));... % 1
                (1+x(1,:)).*(1-x(2,:)).*(1-x(3,:));...
                (1-x(1,:)).*(1+x(2,:)).*(1-x(3,:));... % 3
                (1+x(1,:)).*(1+x(2,:)).*(1-x(3,:));...
                (1-x(1,:)).*(1-x(2,:)).*(1+x(3,:));... % 5
                (1+x(1,:)).*(1-x(2,:)).*(1+x(3,:));...
                (1-x(1,:)).*(1+x(2,:)).*(1+x(3,:));... % 7
                (1+x(1,:)).*(1+x(2,:)).*(1+x(3,:));]/8;
        end
        
        function dNx = gradN(~, x)
            dNx = [-(1-x(2,:)).*(1-x(3,:)) -(1-x(1,:)).*(1-x(3,:)) -(1-x(1,:)).*(1-x(2,:));... % 1
                (1-x(2,:)).*(1-x(3,:)) -(1+x(1,:)).*(1-x(3,:)) -(1+x(1,:)).*(1-x(2,:));...
                -(1+x(2,:)).*(1-x(3,:)) (1-x(1,:)).*(1-x(3,:)) -(1-x(1,:)).*(1+x(2,:));... % 3
                (1+x(2,:)).*(1-x(3,:))  (1+x(1,:)).*(1-x(3,:)) -(1+x(1,:)).*(1+x(2,:));...
                -(1-x(2,:)).*(1+x(3,:)) -(1-x(1,:)).*(1+x(3,:)) (1-x(1,:)).*(1-x(2,:));... % 5
                (1-x(2,:)).*(1+x(3,:)) -(1+x(1,:)).*(1+x(3,:)) (1+x(1,:)).*(1-x(2,:));...
                -(1+x(2,:)).*(1+x(3,:)) (1-x(1,:)).*(1+x(3,:)) (1-x(1,:)).*(1+x(2,:));... % 7
                (1+x(2,:)).*(1+x(3,:))  (1+x(1,:)).*(1+x(3,:)) (1+x(1,:)).*(1+x(2,:))]/8;
        end
    end
    
    methods(Static)
        function res = test_TrilinearBasisFun
            q = fem.HexahedronTrilinear;
            res = fem.BaseFEM.test_BasisFun(q);
            
            % test for correct basis function values on nodes
            [X,Y,Z] = ndgrid(-1:2:1,-1:2:1,-1:2:1);
            p = [X(:) Y(:) Z(:)]';
            res = res && isequal(q.N(p),eye(8));
        end
    end
    
end

