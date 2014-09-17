classdef Model < muscle.Model
% Model: 
%
% Features:
% - regionen mit summe der gewichtungen kleiner 1 ist region mit fett?!
% Scenarios:
% - breiter block an einer seite fest mit allen fasertypen über die länge drin & entsprechendem
% muskelverhalten
% - block mit links fest, rechts per neumann dran ziehen, dann
% spindel-feedback bei zu starker dehnung sorgt für konktraktion
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
                conf = fullmuscle.DebugConfig;
            end
            this = this@muscle.Model(conf);
            
            this.System = fullmuscle.System(this);
            
            % Call the config-specific model configuration
            conf.Model = this;
            conf.configureModel(this);
        
            this.DefaultMu = [1 0 0 3]';
            this.DefaultInput = 1;
            
            % Set the config to the model, triggering geometry related
            % pre-computations
            this.setConfig(conf);
        end
    end
    
end