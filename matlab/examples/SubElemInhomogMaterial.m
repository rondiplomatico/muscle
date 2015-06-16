classdef SubElemInhomogMaterial < muscle.AModelConfig
        
    methods
        function this = SubElemInhomogMaterial(varargin)
            % Creates a Debug simple muscle model configuration.
            this = this@muscle.AModelConfig(varargin{:});
            this.addOption('Flip',1);
            this.init;
            
            this.VelocityBCTimeFun = tools.ConstantUntil(20);
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            
            m.T = 40;
            m.dt = .5;
            m.ODESolver.RelTol = 1e-6;
            m.ODESolver.AbsTol = 1e-3;
            m.DefaultMu(2) = 0;
        end
        
        function tmr = getTendonMuscleRatio(this, points)
            % Returns the [0,1] ratio between tendon and muscle at all
            % specified points
            %
            % This method simply returns an all-zero ratio, meaning muscle only. 
            tmr = zeros(1,size(points,2));
%             tmr(:)=1;
            switch this.Options.Flip
                case 3
                    tmr(points(3,:)<.5 & points(1,:)>.5) = 1;
                    tmr(points(3,:)>2.5 & points(1,:)<.5) = 1;
                case 2
                    tmr(points(2,:)<.5 & points(3,:)>.5) = 1;
                    tmr(points(2,:)>2.5 & points(3,:)<.5) = 1;
                otherwise
                    tmr(points(1,:)<.5 & points(3,:)>.5) = 1;
                    tmr(points(1,:)>2.5 & points(3,:)<.5) = 1;
            end
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            % Single cube with same config as reference element
            xr = linspace(0,3,this.Options.GeoNr+1);
            switch this.Options.Flip
                case 3
                    [pts, cubes] = geometry.Cube8Node.DemoGrid([0 1],[0 1], xr);
                case 2
                    [pts, cubes] = geometry.Cube8Node.DemoGrid([0 1],xr,[0 1]);
                otherwise
                    [pts, cubes] = geometry.Cube8Node.DemoGrid(xr,[0 1],[0 1]);
            end
            %pts(1,:) = pts(1,:)-.5;
            geo = geometry.Cube8Node(pts, cubes);
            geo = geo.toCube27Node;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.PosFE.Geometry;
            switch this.Options.Flip
                case 3
                    displ_dir(:,geo.Elements(1,geo.MasterFaces(5,:))) = true;
                    %displ_dir([2 3],geo.Elements(this.Options.GeoNr,geo.MasterFaces(6,:))) = true;
                case 2
                    displ_dir(:,geo.Elements(1,geo.MasterFaces(3,:))) = true;
                    %displ_dir([2 3],geo.Elements(this.Options.GeoNr,geo.MasterFaces(4,:))) = true;
                otherwise
                    displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
                    %displ_dir([2 3],geo.Elements(this.Options.GeoNr,geo.MasterFaces(2,:))) = true;
            end
            
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            % Determines the dirichlet velocities.
            %
            % The unit for the applied quantities is [mm/ms] = [m/s]
            geo = this.PosFE.Geometry;
            switch this.Options.Flip
                case 3
                    velo_dir(3,geo.Elements(this.Options.GeoNr,geo.MasterFaces(6,:))) = true;
                case 2
                    velo_dir(2,geo.Elements(this.Options.GeoNr,geo.MasterFaces(4,:))) = true;
                otherwise
                    velo_dir(1,geo.Elements(this.Options.GeoNr,geo.MasterFaces(2,:))) = true;
            end
            velo_dir_val(velo_dir) = .02;
        end
        
        function anull = seta0(this, anull)
            anull(this.Options.Flip,:,:) = 1;
        end
    end
    
    methods(Static)
        function res = test_SubElem(gn)
            if nargin < 1
                gn = 3;
            end
            %% 
            v = zeros(3,1);
            o = ones(3,1);
            for fl = 1:3
                c = SubElemInhomogMaterial('GeoNr',gn,'Flip',fl);
                m = c.createModel;
                %m.setGaussIntegrationRule(6);
                [t,y] = m.simulate;
                df = m.getResidualForces(t,y);
                idx = m.getPositionDirichletBCFaceIdx(1,1+(fl-1)*2);
                force = sum(df(idx,end));
                fprintf('Flip %d - Forces at pos dirichlet side at end: %g\n',fl,force);
                v(fl) = force;
                m.plot(t,y,'DF',df,'F',2);
            end
            err = tril(v*o'-(v*o')')
            rel1 = err ./ (v*o')
            rel2 = err ./ (v*o')'
            res = all(abs(rel1(:)) < 1e-4) && all(abs(rel2(:)) < 1e-4)
        end
    end
    
end

