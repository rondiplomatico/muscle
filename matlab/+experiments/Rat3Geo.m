classdef Rat3Geo < handle
    % Contains methods and information about the rat3 geometry
    
    properties(Access=private)
        rady;
        radz;
        ep = .2; %[mm]
        inner = .6;
        lowerflatfactor = .7;
        fibredata;
    end
    
    methods
        
        function this = Rat3Geo
            cent = 1;
            k = kernels.GaussKernel(11.5);
            k2 = kernels.GaussKernel(12.5);
            this.rady = @(x)k.evaluate(x,cent)*4.1;
            this.radz = @(x)k2.evaluate(x,cent).^2*3;
        end
        
        function load(this)
            base = fileparts(which(mfilename('fullpath')));
            this.fibredata = load(fullfile(base,'fibredirections.mat'));
        end
        
        function geo = getGeometry(this)
            % Creates a Belly-shaped muscle geometry with a thin outer
            % layer for apneurosis resolution.
            %
            % The geometry is (by eye-ball!!) fitted to the rat3 geometry
            % from A. Tomalka
            rad = @(x)[this.rady(x) this.radz(x)]';
            geo = Belly.getBelly(5,[-17 25],'Radius',rad,'Layers',[.8 1],'InnerRadius',this.inner);
            % Flatten the lower part
            bot = geo.Nodes(3,:);
            bot(geo.Nodes(3,:) < 0) = bot(geo.Nodes(3,:) < 0)*this.lowerflatfactor;
            % Swap x-y axis
            yx = [geo.Nodes(2,:); geo.Nodes(1,:); bot];
            geo = geometry.Cube27Node(yx,geo.Elements);
        end
        
        function tmr = getTMR(this, x, y, z)
            % Computes the tendon-muscle ratio.
            %
            % This script makes use of the analytical construction of the
            % ellipsoid muscle belly using two separate axis length (+
            % inner radius).
            % The tmr-function is computed based on the distance of any
            % point (x,y,z) to the surface of the belly, as the apneurosis
            % is on the surface.
            %
            % The main axis orientation is x, and y,z are the cross-section
            % axis of the muscle body.
            dist = (y./(this.rady(x)+this.inner)').^2 + (z./(this.radz(x)+this.inner)').^2 - 1;
            dist = dist/this.ep+1;
            upper = dist;
            % Directly set the ratio in some areas
            upper(x >= -2 | z < -.5 | abs(y) > 3) = 0;
            upper(x < -15) = 1;
            
            % Have to re-compute the distance as the geometrie's lower part
            % has been squeezed by the lowerflatfactor.
            dist = (y./(this.rady(x)+this.inner)').^2 + ...
                (z./(this.lowerflatfactor*(this.radz(x)+this.inner))').^2 - 1;
            dist = dist/this.ep+1;
            % Directly set the ratio in some areas
            lower = dist;
            lower(x <= 4 | z > .5 | abs(y) > 3 ) = 0;
            lower(x > 18) = 1;
            tmr = max(min(upper + lower,1),0);            
        end
        
        function a0 = getA0(this, points)
            if isempty(this.fibredata)
                error('no fibre directions loaded. please call "load"');
            end
            pts = this.fibredata.points;
            d = this.fibredata.directions;
            % Swap x/y coordinates here as the raw geometry data is along the y-axis,
            % but the geomery is oriented along the x axis here
            si = scatteredInterpolant(pts,d(:,2));
            dirx = si(points');
            si = scatteredInterpolant(pts,d(:,1));
            diry = si(points');
            si = scatteredInterpolant(pts,d(:,3));
            dirz = si(points');
            a0 = [dirx diry dirz]';
            % Set values at tendon region manually - no fibre direction
            % data there
            manual = points(1,:) < -15 | points(1,:) > 18;
            a0(1,manual) = 0;
            a0(2,manual) = 1;
            a0(3,manual) = 0;
        end
    end
    
    methods(Static)
        function createPlots(dir)
            c = experiments.IsometricActivation('GeoNr',4);
            if nargin < 1
                dir = c.OutputDir;
            end
            m = c.createModel;
            pm = PlotManager;
            pm.UseFileTypeFolders = false;
            pm.ExportDPI = 200;
            pm.SaveFormats = {'jpg','png'};
            pm.NoTitlesOnSave = true;
            m.plotGeometrySetup('PM',pm);
            view([0 0]);
            pm.FilePrefix = 'XZ';
            pm.savePlots(dir);
            
            view([0 90]);
            pm.FilePrefix = 'XY';
            pm.savePlots(dir);
            
            view([90 0]);
            pm.FilePrefix = 'YZ';
            pm.savePlots(dir);
            
            view([38 38]);
            pm.FilePrefix = '';
            pm.savePlots(dir);
        end
    end
    
end

