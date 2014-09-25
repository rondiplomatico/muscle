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
% - dynamische aktivierung von aussen
% - dynamische aktivierung von aussen mit dynamischem zug
% - zyklischer zug
%
% Fragen:
% Sollte der afferent-faktor zur signalübertragung für alle fasertypen gleich sein?
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
    
    methods
        function this = Model(conf)
            if nargin < 1
                conf = fullmuscle.DebugConfig;
            end
            this = this@muscle.Model(conf);
            
            this.System = fullmuscle.System(this);
            
            this.DefaultMu = [1 0 0 3]';
            this.DefaultInput = 1;
            
            % Der tut auch wunderbar :-)
            slv = solvers.MLWrapper(@ode15s);
%             slv.odeopts = odeset('RelTol',1e-9,'AbsTol',1e-9);
            
%             slv = solvers.MLode15i(this);
%             slv.RelTol = .001;
%             slv.AbsTol = .01;
            this.ODESolver = slv;
            
            % Call the config-specific model configuration
            conf.Model = this;
            conf.configureModel(this);
        
            % Set the config to the model, triggering geometry related
            % pre-computations
            this.setConfig(conf);
            
            this.System.prepareSimulation(this.DefaultMu, this.DefaultInput);
        end
    end
    
end