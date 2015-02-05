classdef AExperimentModelConfig < muscle.AModelConfig
    %AEXPERIMENTMODELCONFIG Model configurations for a set of experiments
    
    properties
        OutputDir = []; 
        ImgDir;
        
        % Flag that indicates
%         ICCompMode;
    end
    
    properties(Dependent)
        CurrentConfigNr;
    end
    
    properties(SetAccess=protected)
        NumConfigurations;
        NumOutputs;
        
        % The experimentally determined output values.
        % Must be a NumConfigurations x NumOutputs vector, if set.
        TargetOutputValues;
        
%         HasICComputation = false;
    end
    
    properties(Access=private)
        fCurConfNr = 1;
    end
    
    methods
        function this = AExperimentModelConfig(varargin)
            % Override in subclasses and set NumConfigurations to the
            % number of possible experiment runs with different IC/BCs.
            this = this@muscle.AModelConfig(varargin{:});
        end
        
        function set.CurrentConfigNr(this, nr)
            % Sets the configuration number.
            %
            % Use this in every overridden method to further specify
            % different behaviour
            if nr < 1 || nr > this.NumConfigurations
                error('Please choose one of the %d possible configurations.',this.NumConfigurations);
            end
            this.fCurConfNr = nr;
        end
        
        function value = get.CurrentConfigNr(this)
            value = this.fCurConfNr;
        end
        
%         function x0 = getX0(this, x0)
%             if ~this.ICCompMode
%                 s = load(fullfile(QuickRelease.OutputDir,sprintf('geo%d_ic.mat',this.GeoNr)));
%                 x0 = s.x0;
%             end
%         end
    end
    
    methods(Access=protected)
        function init(this)
            init@muscle.AModelConfig(this);
            
            % Also init directories to reasonable defaults
            if isempty(this.OutputDir)
                mc = metaclass(this);
                [p,n] = fileparts(which(mc.Name));
                outdir = fullfile(p,n);
                this.OutputDir = outdir;
            end
        end
    end
    
    methods
        function set.OutputDir(this, value)
            Utils.ensureDir(value);
            this.OutputDir = value;
            this.ImgDir = fullfile(value,'img');%#ok
        end
    end
    
    methods(Abstract)
        o = getOutputOfInterest(this, t, y);
    end
    
end

