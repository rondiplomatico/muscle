classdef ExperimentRunner < handle    
    
    properties
        RunParallel;
    end
    
    properties(Dependent)
        StoreExperimentData;
    end
    
    properties(SetAccess=private)
        Config;
        Model;
    end
    
    properties(Access=private)
        multipleexperiments = false;
        cache;
        cachedata;
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
            this.StoreExperimentData = true;
        end
        
        function [out, y, ct] = getCachedExperimentConfig(this, nr, mu)
            if nargin < 3
                mu = this.Model.DefaultMu;
            end
            % The cache key is param + config nr
            key = [mu; nr];
            out = [];
            y = [];
            ct = [];
            if this.cachedata && this.cache.hasData(key);
                data = this.cache.getData(key);
                out = data.out;
                y = data.y;
                ct = data.ct;
            end
        end
        
        function [out, y, ct] = runExperimentConfig(this, nr, mu)
            % Runs a single experiment configuration.
            % Coded here for convenience only.
            
            m = this.Model;
            if nargin < 3
                mu = m.DefaultMu;
            end
            [out, y, ct] = this.getCachedExperimentConfig(nr, mu);
            if isempty(out)
                c = this.Config;
                out = NaN(1,c.NumOutputs);
                try
                    % Apply configuration set within ModelConfig
                    % Updates the associated model automatically (precomputations)
                    c.CurrentConfigNr = nr;
                    
                    % Run simulation
                    [t,y,ct] = m.simulate(mu);
                    
                    % Get output of interest
                    tic;
                    out = c.getOutputOfInterest(t,y);
                    ct = ct + toc;
                    
                    if this.cachedata
                        % The cache key is param + config nr
                        key = [mu; nr];
                        data = struct('y',y,'t',t,'mu',mu,'nr',nr,'ct',ct,'out',out);
                        this.cache.addData(key,data);
                    end
                catch ME
                    % Set to NaN to notify there's something wrong
                    out(:) = NaN;
                    ME.getReport
                end
            end
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
                fprintf('Running %d experiment configurations in parallel...\n',c.NumConfigurations);
                closeafterrun = false;
                if matlabpool('size') == 0
                    matlabpool open;
                    closeafterrun = true;
                end
                parfor nr = 1:c.NumConfigurations
                    t = getCurrentTask();
                    fprintf('Worker %d: Running configuration %d\n',t.ID,nr);
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
            this.multipleexperiments = nex > 1;
            
            if this.RunParallel && nex > 1
                closeafterrun = false;
                if matlabpool('size') == 0
                    matlabpool open;
                    closeafterrun = true;
                end
                fprintf('Running %d experiments (%d configs each) in parallel...\n',nex,c.NumConfigurations);
                parfor k=1:nex
                    t = getCurrentTask();
                    fprintf('Worker %d: Running experiment %d\n',t.ID,k);
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
                e = this;
                if nargin < 2
                    [o, ct] = this.runExperiment;%#ok
                    save(fi,'m','o','c','ct','e');
                else
                    [o, ct] = this.runExperiments(mus);%#ok
                    save(fi,'m','o','c','mus','ct','e');
                end
            end
            res = load(fi);
        end
        
        function set.StoreExperimentData(this, value)
            if value && isempty(this.cache)
                dir = fullfile(this.Config.OutputDir,'trajectories',this.Config.getOptionStr(false));
                this.cache = data.FileDataCollection(dir);
            elseif ~isempty(value) && ~value
                this.cache = [];
            end
            this.cachedata = value;
        end
        
        function value = get.StoreExperimentData(this)
            value = this.cachedata;
        end
    end
    
end

