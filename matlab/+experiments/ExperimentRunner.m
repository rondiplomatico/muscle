classdef ExperimentRunner < handle
    %EXPERIMENTRUNNER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        RunParallel;
    end
    
    properties(SetAccess=private)
        Config;
        Model;
    end
    
    properties(Access=private, Transient)
        multipleexperiments = false;
    end
    
    methods
        function this = ExperimentRunner(model)
            this.Model = model;
            config = model.Config;
            if ~isa(config, 'experiments.AExperimentModelConfig')
                error('"config" parameter must be a AExperimentModelConfig');
            end
            this.Config = config;
        end
        
        function out = runExperimentConfig(this, nr, mu)
            % Runs a single experiment configuration.
            % Coded here for convenience only.
            m = this.Model;
            if nargin < 3
                mu = m.DefaultMu;
            end
            c = this.Config;
            % Apply configuration set within ModelConfig
            c.setConfiguration(nr);
            % Update the model (precomputations)
            m.setConfig(c);
            % Run simulation
            out = NaN(1,c.NumOutputs);
%             if this.RunParallel
                try
                    [t,y] = m.simulate(mu);
                    % Get output of interest
                    out = c.getOutputOfInterest(t,y);
                catch ME
                    % Set to NaN to notify there's something wrong
                    out(:) = NaN;
                    ME.getReport
                end
%             else
%                 [t,y] = m.simulate(mu);
%                 % Get output of interest
%                 out = c.getOutputOfInterest(t,y);
%             end
        end
        
        function out = runExperiment(this, mu)
            % Runs the complete set of experiments for a given parameter mu
            if nargin < 2
                mu = this.Model.DefaultMu;
            end
            c = this.Config;
            out = zeros(c.NumConfigurations,c.NumOutputs);
            % Run all the settings!
            if ~this.multipleexperiments
                pi = ProcessIndicator('Running experiment configurations',c.NumConfigurations);
            end
            for nr = 1:c.NumConfigurations
                out(nr,:) = this.runExperimentConfig(nr, mu);
                if ~this.multipleexperiments
                    pi.step;
                end
            end
            if ~this.multipleexperiments
                pi.stop;
            end
        end
        
        function allout = runExperiments(this, mus)
            % Runs the experiments in parallel for all given mu values
            c = this.Config;
            nex = size(mus,2);
            allout = zeros(c.NumConfigurations,c.NumOutputs,nex);
            this.multipleexperiments = true;
            
            if this.RunParallel
                closeafterrun = false;
                if matlabpool('size') == 0
                    matlabpool open;
                    closeafterrun = true;
                end
                parfor k=1:nex
                    allout(:,:,k) = this.runExperiment(mus(:,k));%#ok
                end
                if closeafterrun
                    matlabpool close;
                end
            else
                pi = ProcessIndicator('Running experiment configurations (%d each) for %d parameters',...
                    nex,false,c.NumConfigurations,nex);
                for k=1:nex
                    allout(:,:,k) = this.runExperiment(mus(:,k));
                    pi.step;
                end
                pi.stop;
            end
            % Put number of experiment to first row
            allout = permute(allout,[3 1 2]);
            
            this.multipleexperiments = false;
        end
    end
    
end

