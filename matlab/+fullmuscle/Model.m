classdef Model < muscle.Model
% Model: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-09-16
%
% @new{0,7,dw,2014-09-16} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % The different discrete fibre types this full muscle knows
        FibreTypes = 0:.2:1;
    end
    
    methods
        function this = Model(conf)
            if nargin < 1
                conf = muscle.DebugConfig;
            end
            this = this@muscle.Model(conf);
            
            this.System = fullmuscle.System(this);
            
            % Call the config-specific model configuration
            conf.Model = this;
            conf.configureModel(this);
            
            % Set the config to the model, triggering geometry related
            % pre-computations
            this.setConfig(conf);
        end
    end
    
end