classdef AExperimentModelConfig < muscle.AModelConfig
    %AEXPERIMENTMODELCONFIG Model configurations for a set of experiments
    
    properties
        % Flag that indicates
%         ICCompMode;
    end
    
    properties(SetAccess=protected)
        NumConfigurations;
        NumOutputs;
        
        % The experimentally determined output values.
        % Must be a NumConfigurations x NumOutputs vector, if set.
        TargetOutputValues;
        
%         HasICComputation = false;
    end
    
    properties(SetAccess=private, GetAccess=protected)
        CurrentConfigNr = 1;
    end
    
    methods
        function this = AExperimentModelConfig(geo)
            % Override in subclasses and set NumConfigurations to the
            % number of possible experiment runs with different IC/BCs.
            if nargin < 1
                %[pts, cubes] = geometry.Cube8Node.DemoGrid;
                geo = geometry.Cube27Node;
            end
            
            this = this@muscle.AModelConfig(geo);
        end
        
        function setConfiguration(this, nr)
            % Sets the configuration number.
            %
            % Use this in every overridden method to further specify
            % different behaviour
            this.CurrentConfigNr = nr;
        end
        
%         function x0 = getX0(this, x0)
%             if ~this.ICCompMode
%                 s = load(fullfile(QuickRelease.OutputDir,sprintf('geo%d_ic.mat',this.GeoNr)));
%                 x0 = s.x0;
%             end
%         end
    end
    
    methods(Abstract)
        o = getOutputOfInterest(this, t, y);
    end
    
end

