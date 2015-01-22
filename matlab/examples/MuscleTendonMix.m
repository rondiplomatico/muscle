classdef MuscleTendonMix < muscle.AModelConfig
    % Muscle - Tendon mixed geometries
    %
    % Contains several examples:
    % Variant 1: A 4 element long block with 4 block-wise constant,
    % linearly changing muscle-to-tendon material. The tendon part is at a
    % dirichlet xz-fixed side and the other side gets increasingly pulled over
    % 4 sec until 4000kPa
    %
    % Variant 2: A single element with the top 9 gauss points as tendons
    % (tmr=1) and the lower part muscle (tmr=0). Fixed (xz only) on right
    % side and pulled on left side over 4 seconds increasing to 4000kPa
    %
    % Variant 3: A 2x2x3 geometry on domain [0 2]x[0 6]x[0 2] with a
    % nonlinear muscle/tendon function (use plotTMRFun to see). Fixed xz on
    % right side (with xyz fix in center node) and pulled on left with
    % 4000kPa. The fibres are at 45deg throughout the muscle
    %
    % Variant 4: The same geometry as Variant 4 with same muscle/tendon
    % ratio, but this time completely fixed at both long ends. An
    % activation ramp is applied over 4s from 0 to 1.
    %
    % Variant 5: The same geometry as variant 2 (single cube) with same
    % setup. No pulling force is applied but the muscle activated over 4s
    % to 1
    %
    % Variant 6: A belly-shaped 
    
    properties
        Variant;
    end
    
    properties(SetAccess=private)
        ylen;
    end
    
    methods
        function this = MuscleTendonMix(variantnr)
            if nargin < 1
                variantnr = 1;
            end
            switch variantnr
                case 1
                    % Four cubes in a row
                    [pts, cubes] = geometry.Cube8Node.DemoGrid(0:1,0:4,0:1);
                case {2 5}
                    % Single cube with same config as reference element
                    [pts, cubes] = geometry.Cube8Node.DemoGrid(0:1,0:1,0:1);
                case {3 4}
                    % 2x4x2 geometry
                    [pts, cubes] = geometry.Cube8Node.DemoGrid(0:2,0:2:6,0:2);
            end
            if any(variantnr == [6 7])
                yl = 10;
                geo = Belly.getBelly(6,yl,1,.2,2);
            else
                geo = geometry.Cube8Node(pts, cubes);
                geo = geo.toCube27Node;
            end
            this = this@muscle.AModelConfig(geo);
            this.Variant = variantnr;
            if any(variantnr == [6 7])
                this.ylen = yl;
            end
            %this.NeumannCoordinateSystem = 'global';
        end
        
        function configureModel(this, m)
            m.T = 4;
            m.dt = .01;
            m.DefaultInput = 1;
            f = m.System.f;
            f.alpha = @(t)0;
            f.Pmax = 250;
            
            os = m.ODESolver;
            os.RelTol = .0001;
            os.AbsTol = .05;
            
            mu = [1; 0; 0; 0
                    %% Anisotropic parameters muscle+tendon (Markert)
                    4.02; 38.5; 7990; 16.6
                    %% Isotropic parameters muscle+tendon (Moonley-Rivlin)
                    35.6; 3.86; 2310; 1.15e-3]; % Micha
                    % 6.352e-10; 3.627 [Kpa] thomas alt
            switch this.Variant
                case {1 2}
                    mu(3) = 4000;
                case 3
                    mu(3) = 4000;
                    m.ODESolver.AbsTol = .1;
                case 4
                    mu(2) = m.T;
                    mu(3) = 0;
                    m.ODESolver.RelTol = .1;
                    m.ODESolver.AbsTol = .1;
                case 5
                    mu(2) = m.T;
                    m.ODESolver.RelTol = .1;
                    m.ODESolver.AbsTol = .1;
                case 6
                    mu(3) = 10000;
                case 7
                    m.T = 30;
                    m.dt = .05;
                    mu(2) = m.T*2/3; % activate 2/3rds of the time
                    mu(3) = 0;
            end
            m.DefaultMu = mu;
            
            % Ramp up the external pressure
            m.System.Inputs{1} = this.getAlphaRamp(1,1);
        end
        
        function tmr = getTendonMuscleRatio(this, points)
            % Returns the [0,1] ratio between tendon and muscle at all
            % specified points
            %
            % This method simply returns an all-zero ratio, meaning muscle only. 
            tmr = zeros(1,size(points,2));
            switch this.Variant
                case 1
                    ratios = linspace(0,1,4);
                    for k=1:4
                        tmr(points(2,:) >= k-1 .* points(2,:) < k) = ratios(k);
                    end
                case {2 5}
                    tmr(points(3,:) > .66) = 1;
                case {3 4}
                    tmr = this.TMRFun(points(2,:),points(3,:));
                    % Take only half - too stiff otherwise
                    tmr = tmr/2;
                case {6 7}
                    % Set the middle 20% percent to 100% muscle
                    zeroperc = .2;
                    
                    y = points(2,:);
                    cent = this.ylen/2;
                    left = y <= cent*(1-zeroperc/2);
                    tmr(left) = 1-y(left)/(cent*(1-zeroperc/2));
                    right = y > cent + this.ylen*zeroperc/2;
                    tmr(right) = (y(right)-(this.ylen*(.5+zeroperc/2)))/...
                        (this.ylen*(.5-zeroperc/2));
            end
            % Use max .2 tendon - works but needs to be checked for higher
            % tendon parts (and current default parameter mu!!!!)
            if this.Variant == 4
                tmr = 2*tmr/5;
            end
        end
        
        function tmr = TMRFun(~, z, y)
            f = @(z,y)min(1,max(0,-.5+.5./((z/2+.1).*(1.5*y+.08*z+.1))));
            tmr = f(z,y) + f(6-z,2-y);
        end
        
        function plotTMRFun(this)
            [Y,Z] = meshgrid(0:.01:2.1,0:.01:6.1);
            tmr = this.TMRFun(Z,Y);
            pm = PlotManager;
            pm.LeaveOpen = true;
            ax = pm.nextPlot('tmr_func','Muscle/Tendon ratio function');
            surfc(Z,Y,tmr,'Parent',ax);
            axis(ax,'equal');
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            switch this.Variant
                % Pull on front face
                case {1 2 3}
                    if any(elemidx == [1 2 7 8]) && faceidx == 3
                        P = 1;
                    end
                case 6
                    if any(elemidx == 1:4) && faceidx == 3
                        P = 1;
                    end
            end
        end
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            
            switch this.Variant
                case {1 2 5}
                    % Fix all on left and only the y,z directions of the back right
                    % nodes
                    displ_dir(2,geo.Elements(geo.NumElements,geo.MasterFaces(4,:))) = true;
                case 3
                    displ_dir(2,geo.Elements([5 6 11 12],geo.MasterFaces(4,:))) = true;
                    displ_dir(:,geo.Elements(5,geo.MasterFaces(4,9))) = true;
                case 4
                    %displ_dir(2,geo.Elements([5 6 11 12],geo.MasterFaces(4,:))) = true;
                    %displ_dir(2,geo.Elements([1 2 7 8],geo.MasterFaces(3,:))) = true;
                    displ_dir(:,geo.Elements([5 6 11 12],geo.MasterFaces(4,:))) = true;
                    displ_dir(:,geo.Elements([1 2 7 8],geo.MasterFaces(3,:))) = true;
                case 6
                    displ_dir(:,geo.Elements(13:16,geo.MasterFaces(4,:))) = true;
                case 7
                    displ_dir(:,geo.Elements(1:4,geo.MasterFaces(3,:))) = true;
                    last = size(geo.Elements,1)-3:size(geo.Elements,1);
                    displ_dir(:,geo.Elements(last,geo.MasterFaces(4,:))) = true;
            end
        end
        
        function anull = seta0(this, anull)
            switch this.Variant
                case {1 2 5 6 7}
                    % Direction is y
                    anull(2,:,:) = 1;
                case {3 4}
                    %anull(2,:,:) = 1;
                    anull([2 3],:,:) = -1;
            end
        end
    end
    
end

