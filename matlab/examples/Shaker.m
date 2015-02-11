classdef Shaker < muscle.AModelConfig
    
    properties
        ylen;
    end
    
    methods
        function this = Shaker(varargin)
            this = this@muscle.AModelConfig(varargin{:});
            this.init;
            this.VelocityBCTimeFun = tools.Sinus(50); % 50Hz
        end
        
        function configureModel(this, m)
            m.T = 40;
            m.dt = .05;
            
            mu = m.DefaultMu;
            % Small viscosity
            mu(1) = 1e-5;
            m.DefaultMu = mu;
        end
        
        function tmr = getTendonMuscleRatio(this, points)
            % Returns the [0,1] ratio between tendon and muscle at all
            % specified points
            %
            % This method simply returns an all-zero ratio, meaning muscle only. 
            tmr = zeros(1,size(points,2));
            
            % Set the middle 20% percent to 100% muscle
            zeroperc = .2;

            y = points(2,:);
            cent = this.ylen/2;
            left = y <= cent*(1-zeroperc/2);
            tmr(left) = 1-y(left)/(cent*(1-zeroperc/2));
            right = y > cent + this.ylen*zeroperc/2;
            tmr(right) = (y(right)-(this.ylen*(.5+zeroperc/2)))/...
                (this.ylen*(.5-zeroperc/2));
            
            tmr = 2*tmr/5;
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            belly = Belly.getBelly(4,10,'Radius',1,'InnerRadius',.1,'Gamma',2);
            geo = belly.scale(10);
            this.ylen = 100;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            % Fix ends in xz direction
            displ_dir([1 3],geo.Elements(1:4,geo.MasterFaces(3,:))) = true;
            displ_dir([1 3],geo.Elements(13:16,geo.MasterFaces(4,:))) = true;
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            % Fix ends in xz direction
            velo_dir(2,geo.Elements(1:4,geo.MasterFaces(3,:))) = true;
            velo_dir(2,geo.Elements(13:16,geo.MasterFaces(4,:))) = true;
%             for k = [1 2 7 8 5 6 11 12]  
%                 pos = geo.Elements(k,geo.MasterFaces(3,:));
%                 velo_dir(1,pos) = true;
%                 %velo_dir_val(1,pos) = 1;
%             end
            velo_dir_val(velo_dir) = 1;
        end
        
        function anull = seta0(~, anull)
            % Direction is y
            anull(2,:,:) = 1;
        end
    end
    
    methods(Static)
        function test_Shaker
            m = muscle.Model(Shaker);
            m.simulateAndPlot;
        end
    end
    
end

