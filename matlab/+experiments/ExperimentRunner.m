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
    
    properties(Access=private)
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
            
            % Do parallel if parallel toolbox is around!
            this.RunParallel = exist('matlabpool','file') == 2;
        end
        
        function [out, y, ct] = runExperimentConfig(this, nr, mu)
            % Runs a single experiment configuration.
            % Coded here for convenience only.
            
            m = this.Model;
            if nargin < 3
                mu = m.DefaultMu;
            end
            c = this.Config;
            % Apply configuration set within ModelConfig
            c.CurrentConfigNr = nr;
            % Update the model (precomputations)
            m.setConfig(c);
            % Run simulation
            out = NaN(1,c.NumOutputs);
%             if this.RunParallel
                try
                    [t,y,ct] = m.simulate(mu);
                    % Get output of interest
                    tic;
                    out = c.getOutputOfInterest(t,y);
                    ct = ct + toc;
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
        
        function [out, ct] = runExperiment(this, mu)
            % Runs the complete set of experiments for a given parameter mu
            if nargin < 2
                mu = this.Model.DefaultMu;
            end
            c = this.Config;
            out = zeros(c.NumConfigurations,c.NumOutputs);
            ct = zeros(c.NumConfigurations,1);
            % Run all the settings!
            if this.RunParallel && ~this.multipleexperiments
                fprintf('Running %d experiment configurations in parallel...',c.NumConfigurations);
                if matlabpool('size') == 0
                    matlabpool open;
                    closeafterrun = true;
                end
                parfor nr = 1:c.NumConfigurations
                    [out(nr,:), ~, ct(nr)] = this.runExperimentConfig(nr, mu);%#ok
                end
                if closeafterrun
                    matlabpool close;
                end
                fprintf('done');
            else
                if ~this.multipleexperiments
                    pi = ProcessIndicator('Running experiment configurations',c.NumConfigurations);
                end
                for nr = 1:c.NumConfigurations
                    [out(nr,:), ~, ct(nr)] = this.runExperimentConfig(nr, mu);
                    if ~this.multipleexperiments
                        pi.step;
                    end
                end
                if ~this.multipleexperiments
                    pi.stop;
                end
            end
        end
        
        function [allout, allct] = runExperiments(this, mus)
            % Runs the experiments in parallel for all given mu values
            c = this.Config;
            nex = size(mus,2);
            allout = zeros(c.NumConfigurations,c.NumOutputs,nex);
            allct = zeros(c.NumConfigurations,nex);
            this.multipleexperiments = true;
            
            if this.RunParallel
                closeafterrun = false;
                if matlabpool('size') == 0
                    matlabpool open;
                    closeafterrun = true;
                end
                fprintf('Running %d experiments (%d configs each) in parallel...',nex,c.NumConfigurations);
                parfor k=1:nex
                    [allout(:,:,k), allct(:,k)] = this.runExperiment(mus(:,k));%#ok
                end
                fprintf('done');
                if closeafterrun
                    matlabpool close;
                end
            else
                pi = ProcessIndicator('Running experiment configurations (%d each) for %d parameters',...
                    nex,false,c.NumConfigurations,nex);
                for k=1:nex
                    [allout(:,:,k), allct(:,k)] = this.runExperiment(mus(:,k));
                    pi.step;
                end
                pi.stop;
            end
            % Put number of experiment to first row
            allout = permute(allout,[3 1 2]);
            
            this.multipleexperiments = false;
        end
        
        function res = runExperimentsCached(this, mus)
            c = this.Config;
            fi = fullfile(c.OutputDir,[c.getOptionStr '.mat']);
            if ~(exist(fi,'file') == 2)
                m = this.Model;
                if nargin < 2
                    [o, ct] = this.runExperiment;%#ok
                    save(fi,'m','o','c','ct');
                else
                    [o, ct] = this.runExperiments(mus);%#ok
                    save(fi,'m','o','c','mus','ct');
                end
            end
            res = load(fi);
        end
    end
    
end

