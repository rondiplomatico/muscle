classdef CFibreDemo < fullmuscle.AModelConfig
    
    properties(SetAccess=private)
        Version;
    end
    
    methods
        function this = CFibreDemo
            % Creates a Debug simple muscle model configuration.
            %
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid([0 10],0:10:40,[0 10]);
            geo = geometry.Cube8Node(pts, cubes);
            this = this@fullmuscle.AModelConfig(geo.toCube27Node);
            this.Version = version;
        end
        
        function configureModel(this, m)
            configureModel@fullmuscle.AModelConfig(this, m);
            
            m.T = 2500;
            m.dt = 1;
            
            m.DefaultMu(1:4) = [1; 0; 1; 0];
            m.DefaultMu(13) = 250;
            m.EnableTrajectoryCaching = true;
            s = m.System.Plotter;
            s.GeoView = [72 74];
        end
        
        function u = getInputs(this)
            null = @(t)zeros(size(t));
            u(1,1:3) = {null null null};
            u{2,1} = this.getAlphaRamp(10,1);
            u{2,2} = this.getAlphaRamp(200,1);
            u{2,3} = this.getAlphaRamp(500,1);
        end
        
    end
    
    methods(Access=protected)
        
        function ft = getFibreTypes(~) 
            ft = [0 .2 .4 .6 .8 1];
        end
        
        function sp = getSpindlePos(~)
            % No spindles (feedback) here
            sp = [];
        end
        
        function ftw = getFibreTypeWeights(this)
            % Get pre-initialized all zero weights
            ftw = getFibreTypeWeights@fullmuscle.AModelConfig(this);

            ls = LinearSplitOfOne('y',5,40);
            splitfun = ls.getAllFunctionsParsed;
            
            fe = this.PosFE;
            geo = fe.Geometry;
            for k = 1:geo.NumElements
                gp_ypos = geo.Nodes(2,geo.Elements(k,:)) * fe.Ngp(:,:,k);
                ftw(:,:,k) = splitfun(gp_ypos)';
            end
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.PosFE.Geometry;
            % Always fix back side
            for k = 1:geo.NumElements
                displ_dir(:,geo.Elements(k,geo.MasterFaces(1,:))) = true;
            end
        end
        
        function anull = seta0(~, anull)
            anull(1,:,:) = 1;
        end
    end
    
end

