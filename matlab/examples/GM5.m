classdef GM5 < muscle.AModelConfig
    % This AModelConfig involves experiments around the GM5 muscle,
    % experimentally measured by Tobias Siebert.
    %
    % Muscle mass: 898mg = 0.898g
    % Muscle density: 1.056 g/cm³ = 0.001056g/mm³
    % Muscle volume: 850.37mm³
    %
    % Transversal loading experiment:
    % 
    
    properties
        velo_time = 20; % [ms]
        target_dist = -6; % [mm]
    end
    
    methods
        function this = GM5
            
            l = 42; % muscle length [mm]
            h = 9.9; % muscle height [mm] [orig 9.9 bei 45°]
            b = 22; % muscle outmost width [mm]
            
            % Internal: Set the plate for muscle bottom flattening to -2.5
%             plate_at = -2.5;
            plate_at = 0;
            
            geo = GM5.createGeometry(l,b,h,plate_at);
            this = this@muscle.AModelConfig(geo);
            
            this.a0CoordinateSystem = 'reference';
%             this.a0CoordinateSystem = 'master';
        end
        
        function configureModel(this, m)
            % Overload this method to set model-specific quantities like
            % simulation time etc
            
            m.MuscleDensity = 1.056e-6; % [kg/mm³] = 1.056 [g/cm³]
            
            m.System.ApplyVelocityBCUntil = 20;
            
            m.T = 200;
            m.dt = 1;
            m.DefaultMu = [.1; 0];
            m.System.f.alpha = this.getAlphaRamp(50,1,30);
            os = m.ODESolver;
            os.RelTol = .001;
            os.AbsTol = .08;
        end
    end
    
    methods(Access=protected)
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            % Front face
            for e = 1:4
                displ_dir(:,geo.Elements(e,geo.MasterFaces(3,:))) = true;
            end
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            geo = this.PosFE.Geometry;
            for e = 9:12
                velo_dir(2,geo.Elements(e,geo.MasterFaces(4,:))) = true;
            end
            velo_dir_val(velo_dir) = this.target_dist/this.velo_time;
        end
                
        function anull = seta0(~, anull)
            % The elements are aligned so that the fibres go in y-direction
            anull(2,:,:) = 1;
            anull(3,:,:) = .5;
            %anull(1,:,:) = 1;
        end
    end
    
    methods(Static)
        function geo = createGeometry(l, w, h, plate_at)
            % Creates the geometry instance
            
            % Recompute h so that overall height stays at original value
            h = h + plate_at;
            
            % Compute the 
            z = linspace(0,l,600);
            str = 'log(z+1)+fliplr(.45*sqrt(z))';
            fbase = eval(['@(z)' str ';']);
            fzmin = min(fbase(z));%#ok
            fbase = eval(['@(z)' str '-fzmin;']);
            fzmax = max(fbase(z));%#ok
            fbase = eval(['@(z)(' str '-fzmin)/fzmax;']);
            plot(z,fbase(z));

            fz = @(z)[fbase(z)*w/2; fbase(z)*h/2];
            geo = Belly.getBelly(3,l,fz,[2 1.5],5,plate_at);
        end
    end
    
end

