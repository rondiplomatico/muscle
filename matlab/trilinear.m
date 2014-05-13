classdef trilinear < fembase
    %TRILINEAR Summary of this class goes here
    %   Detailed explanation goes here
    %
    
    methods
        function this = trilinear(geo)
            if nargin < 1
                geo = geometry.Cube8Node;
            end
            
            this = this@fembase(geo);
            
            %this.EdgeIndices = 1:8;
            this.init;
        end
        
        function init(this)
            % Trilinear basis functions
            this.N = @(x)[   (1-x(1,:)).*(1-x(2,:)).*(1-x(3,:));... % 1
                (1+x(1,:)).*(1-x(2,:)).*(1-x(3,:));...
                (1-x(1,:)).*(1+x(2,:)).*(1-x(3,:));... % 3
                (1+x(1,:)).*(1+x(2,:)).*(1-x(3,:));...
                (1-x(1,:)).*(1-x(2,:)).*(1+x(3,:));... % 5
                (1+x(1,:)).*(1-x(2,:)).*(1+x(3,:));...
                (1-x(1,:)).*(1+x(2,:)).*(1+x(3,:));... % 7
                (1+x(1,:)).*(1+x(2,:)).*(1+x(3,:));]/8;
            
            this.gradN = @(x)[-(1-x(2,:)).*(1-x(3,:)) -(1-x(1,:)).*(1-x(3,:)) -(1-x(1,:)).*(1-x(2,:));... % 1
                (1-x(2,:)).*(1-x(3,:)) -(1+x(1,:)).*(1-x(3,:)) -(1+x(1,:)).*(1-x(2,:));...
                -(1+x(2,:)).*(1-x(3,:)) (1-x(1,:)).*(1-x(3,:)) -(1-x(1,:)).*(1+x(2,:));... % 3
                (1+x(2,:)).*(1-x(3,:))  (1+x(1,:)).*(1-x(3,:)) -(1+x(1,:)).*(1+x(2,:));...
                -(1-x(2,:)).*(1+x(3,:)) -(1-x(1,:)).*(1+x(3,:)) (1-x(1,:)).*(1-x(2,:));... % 5
                (1-x(2,:)).*(1+x(3,:)) -(1+x(1,:)).*(1+x(3,:)) (1+x(1,:)).*(1-x(2,:));...
                -(1+x(2,:)).*(1+x(3,:)) (1-x(1,:)).*(1+x(3,:)) (1-x(1,:)).*(1+x(2,:));... % 7
                (1+x(2,:)).*(1+x(3,:))  (1+x(1,:)).*(1+x(3,:)) (1+x(1,:)).*(1+x(2,:))]/8;
            
            init@fembase(this);
        end
    end
    
    methods(Static)
        function res = test_TrilinearBasisFun
            q = trilinear;
            res = fembase.test_BasisFun(q);
            
            % test for correct basis function values on nodes
            [X,Y,Z] = ndgrid(-1:2:1,-1:2:1,-1:2:1);
            p = [X(:) Y(:) Z(:)]';
            res = res && isequal(q.N(p),eye(8));
        end
    end
    
end

