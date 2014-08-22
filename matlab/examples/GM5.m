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
    
    methods
        function this = GM5
            
            l = 35; % muscle length [mm]
            h = 9.9; % muscle height [mm]
            b = 25; % muscle outmost width [mm]
            
            geo = GM5.createGeometry(l,b,h);
            this = this@muscle.AModelConfig(geo);
        end
        
        function configureModel(this, m)
            % Overload this method to set model-specific quantities like
            % simulation time etc
            
            m.MuscleDensity = 1.056e-6; % [kg/mm³] = 1.056 [g/cm³]
            
            m.T = 200;
            m.dt = 1;
            m.DefaultMu = [.1; 0];
            m.System.f.alpha = this.getAlphaRamp(50,1);
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
            for e = 9:12
                displ_dir(:,geo.Elements(e,geo.MasterFaces(4,:))) = true;
            end
        end
                
        function anull = seta0(~, anull)
            % The elements are aligned so that the fibres go in x-direction
            anull(1,:,:) = 1;
        end
    end
    
    methods(Static)
        function geo = createGeometry(l, w, h)
            % Creates the geometry instance
            
            % Internal: Set the plate for muscle bottom flattening to -2.5
            plate_at = -2.5;
            % Recompute h so that overall height stays at original value
            h = h + plate_at;
            
            % Compute the 
            z = linspace(0,l,600);
            str = 'log(z+1)+fliplr(.3*sqrt(z))';
            fbase = eval(['@(z)' str ';']);
            fzmin = min(fbase(z));%#ok
            fbase = eval(['@(z)' str '-fzmin;']);
            fzmax = max(fbase(z));%#ok
            fbase = eval(['@(z)(' str '-fzmin)/fzmax;']);
%             plot(z,fbase(z));

            fz = @(z)[fbase(z)*w/2; fbase(z)*h/2];
            geo = Belly.getBelly(3,l,fz,[2 1.5],5,-2.5);
        end
    end
    
end

