classdef Sprenger32Elem < muscle.AModelConfig
    
    methods
        function this = Sprenger32Elem
            s = load(fullfile(fileparts(which(mfilename)),'..','CMISS','Sprenger32Elem.mat'));
            this = this@muscle.AModelConfig(s.geo20);
        end
        
        function configureModel(this, model)
            model.T = 1;
            model.dt = .01;
            f = model.System.f;
            f.System.Viscosity = 1;
            
            %% Material configuration from CMISS/3Elem_sprenger.xml
            % malpha_calculation
            f.alpha = @(t)1; % [-]
            
            % c1M = 3.56463903963e-02 MPa
            f.c10 = 35.6463903963; % [kPa]
            
            % c2M = 3.859558659683e-3 MPa
            f.c01 = 3.859558659683; % [kPa]
            
            % c3M = 3.05907e-10 MPa
            f.b1 = 3.05907e-7; % [kPa]
            
            % c4M 
            f.d1 = 47.270456264135881; % [-]
            
            % sigma_max_calculation = 0.3 in MPa
            f.Pmax = 300; % [kPa]
            
            % lambda_ofl_calculation
            f.lambdafopt = 1.4; % [-]
            
            f.ForceLengthFun = @(ratio)(ratio<=1).*exp(-((1-ratio)/.57).^4) + (ratio>1).*exp(-((ratio-1)/.14).^3);
            % The derivative of the force-length function as function handle
            f.ForceLengthFunDeriv = @(ratio)(ratio<=1).*((1/.57)*(((1-ratio)/.57).^3).*exp(-((1-ratio)/.57).^4)) ...
                - (ratio > 1) .* ((1/.14) .* (((ratio-1)/.14).^2) .* exp(-((ratio-1)/.14).^3));
        end
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            % 1:4 one side
            displ_dir(1,geo.Elements(1,[1 4 6 9 13])) = true;
            displ_dir(:,geo.Elements(1,[11 16 18])) = true;
            displ_dir(1,geo.Elements(2,[1 4 6 11 18])) = true;
            displ_dir(:,geo.Elements(2,[9 13 16])) = true;
            displ_dir(1,geo.Elements(3,[1 9 13 16 18])) = true;
            displ_dir(:,geo.Elements(3,[4 6 11])) = true;
            displ_dir(1,geo.Elements(4,[6 11 13 16 18])) = true;
            displ_dir(:,geo.Elements(4,[1 4 9])) = true;
            % 29:32 other side
%             displ_dir(1,geo.Elements(29,[3 5 8 10 15])) = true;
%             displ_dir(:,geo.Elements(29,[12 17 20])) = true;
%             displ_dir(1,geo.Elements(30,[3 5 8 12 20])) = true;
%             displ_dir(:,geo.Elements(30,[10 13 17])) = true;
%             displ_dir(1,geo.Elements(31,[3 10 15 17 20])) = true;
%             displ_dir(:,geo.Elements(31,[5 8 12])) = true;
%             displ_dir(1,geo.Elements(32,[8 12 15 17 20])) = true;
%             displ_dir(:,geo.Elements(32,[3 5 10])) = true;
        end
        
        function anull = seta0(~, anull)
            % Direction is xz
            anull(1,:,:) = 1;
        end
    end
    
end

