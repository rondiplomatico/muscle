classdef System < muscle.System;
% System: 
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
    
    properties(SetAccess=private)
        num_motoneuron_dof;
        num_sarco_dof;
        num_spindle_dof;
        num_all_dof;
        off_moto;
        off_sarco;
        off_spindle;
        
        input_motoneuron_link_idx;
        moto_input_noise_factors;
        sarco_output_idx;
        
        % The upper limit polynomial for maximum mean current dependent on
        % the fibre type.
        %
        % Used at the assembly of B to provide suitable coefficient
        % functions.
        %
        % See also: models.motoneuron.experiments.ParamDomainDetection
        upperlimit_poly;
        
        % The offset for the sarcomere signals at t=0. Used to create
        % alpha(X,0) = 0
        %
        % See also: assembleX0
        sarco_mech_signal_offset
        
        Motoneuron;
        Spindle;
        Sarcomere;
        HasSpindle;
    end
    
    properties(Access=private)
        nfibres;
    end
    
    methods
        function this = System(model)
            this = this@muscle.System(model);
            this.f = fullmuscle.Dynamics(this);
            
            % Load mean current limiting polynomial
            s = load(models.motoneuron.Model.FILE_UPPERLIMITPOLY);
            this.upperlimit_poly = s.upperlimit_poly;
            
            % First row is neumann input
            this.Inputs{1,1} = model.Config.getAlphaRamp(30,1);
            % Second row is external mean current input
%             this.Inputs{2,1} = @(t)1;
            
            this.Motoneuron = models.motoneuron.Motoneuron;
            this.Sarcomere = fullmuscle.Sarcomere;
            
            % Compile information for plotting
            this.Plotter = fullmuscle.MusclePlotter(this);
        end
        
        function configUpdated(this)
            mc = this.Model.Config;
            ft = mc.FibreTypes;
            nf = length(ft);
            this.nfibres = nf;
            
            %% Configure motoneurons
            this.Motoneuron.setType(ft);
            
            %% Configure sarcomeres
            this.Sarcomere.setType(ft);
            
            %% Configure Spindle (if set)
            hassp = ~isempty(mc.SpindlePositions);
            this.HasSpindle = hassp;
            this.Spindle = [];
            if hassp
                this.Spindle = fullmuscle.Spindle;
            end
            
            configUpdated@muscle.System(this);
        end
        
        function setConfig(this, mu, inputidx)
            setConfig@muscle.System(this, mu, inputidx);
            
            if ~isempty(inputidx)
                % Create an input substitute that uses the true external
                % function and creates the effective noisy signal from it
                maxcurrents = polyval(this.upperlimit_poly,this.Model.Config.FibreTypes);
                
                % Get noise from motoneuron class
                no = this.Motoneuron.TypeNoise;
                bno = this.Motoneuron.BaseNoise;
                
                % First row is neumann input
                uneum = this.Inputs{1,inputidx};
                
                % Second row is external mean current input
                uext = this.Inputs{2,inputidx};
                
                ustr = '@(t)[mu(3)*uneum(t); bno(round(t)+1); ';
                for k=1:this.nfibres
                    rowfun = sprintf('no(%d,round(t)+1)*min(%g,mu(4)*uext(t)); ', k, maxcurrents(k));
                    ustr = [ustr rowfun];%#ok
                end
                ustr = [ustr ']'];
                this.u = eval(ustr);
            end
        end
        
        function uvwall = includeDirichletValues(this, t, uvw)
            uvwall_mech = includeDirichletValues@muscle.System(this, t, uvw(1:this.NumTotalDofs,:));
            uvwall = [uvwall_mech; uvw(this.NumTotalDofs+1:end,:)];
        end
        
    end
    
    methods(Access=protected)
        
        function updateDimensions(this, mc)
            updateDimensions@muscle.System(this, mc);
            
            this.num_motoneuron_dof = 6*this.nfibres;
            % Motoneurons are beginning after mechanics
            this.off_moto = this.NumTotalDofs; 

            % Sarcomeres are beginning after motoneurons
            this.num_sarco_dof = 56*this.nfibres;
            this.off_sarco = this.off_moto + this.num_motoneuron_dof;
            
            % Spindles are beginning after sarcomeres
            this.off_spindle = this.off_sarco + this.num_sarco_dof;
            this.num_spindle_dof = 0;
            if this.HasSpindle
                this.num_spindle_dof = 9*this.nfibres;
            end
            this.num_all_dof = this.off_spindle + this.num_spindle_dof;
            
            % Get the positions where the input signal is mapped to the
            % motoneurons
            this.input_motoneuron_link_idx = this.off_moto + (2:6:6*this.nfibres);
            
            this.sarco_output_idx = this.off_sarco + (53:56:56*this.nfibres);
        end
        
        function x0 = assembleX0(this)
            x0 = zeros(this.num_all_dof,1);
            % Get muscle x0
            x0(1:this.NumTotalDofs) = assembleX0@muscle.System(this);
            
            % Load dynamic/affine x0 coefficients for moto/sarco system
            % from file
            mc = metaclass(this);
            s = load(fullfile(fileparts(which(mc.Name)),'x0coeff.mat'));
            x0_motorunit = dscomponents.AffineInitialValue;
            m = size(s.coeff,1);
            for k=1:m
                x0_motorunit.addMatrix(sprintf('polyval([%s],mu(1))',...
                    sprintf('%g ',s.coeff(k,:))),full(sparse(k,1,1,m,1)));
            end
            
            ft = this.Model.Config.FibreTypes;
            smoff = zeros(this.nfibres,1);
            for k=1:this.nfibres
                x0ms = x0_motorunit.evaluate(ft(k));
                % add moto
                x0(this.off_moto + 6*(k-1) + (1:6)) = x0ms(1:6);
                % add sarco
                x0(this.off_sarco + 56*(k-1) + (1:56)) = x0ms(7:end);
                smoff(k) = x0ms(6+53);
                if this.HasSpindle
                    % add spindle
                    x0(this.off_spindle + 9*(k-1) + (1:9)) = this.Spindle.y0;
                end
            end
            this.sarco_mech_signal_offset = smoff;
        end
        
        function Baff = assembleB(this)
            Baff = dscomponents.AffLinInputConv;
            
            %% Add neumann forces B into larger affine matrix in first
            % column
            % Assemble B matrix for neumann conditions from muscle.System
            Baff_neumann = assembleB@muscle.System(this);
            if ~isempty(Baff_neumann)
                B = sparse(this.num_all_dof,this.nfibres+2);
                B(1:this.NumTotalDofs,1) = Baff_neumann.getMatrix(1);
                Baff.addMatrix(Baff_neumann.funStr{1},B);
            end
            
            %% Add motoneuron external input signal
            i = this.input_motoneuron_link_idx;
            s = this.Motoneuron.FibreTypeNoiseFactors;
            % Initialize with base noise entries in second column
            B = sparse(i,ones(size(i))+1,s,this.num_all_dof,this.nfibres+2);
            % Add single contribution for each fibre type thereafter
            for k=1:this.nfibres
                B(i(k),k+2) = s(k);%#ok
            end
            Baff.addMatrix('1',B);
        end
        
        function Daff = assembleDampingMatrix(this)
            Daff_mech = assembleDampingMatrix@muscle.System(this);
            
            Daff = dscomponents.AffLinCoreFun(this);
            extradofs = 6+56;
            if this.HasSpindle
                extradofs = extradofs + 9;
            end
            extra = extradofs*this.nfibres;
            D = blkdiag(Daff_mech.getMatrix(1),sparse(extra,extra));
            Daff.addMatrix('mu(1)',D);
        end
        
        function M = assembleMassMatrix(this)
            M = assembleMassMatrix@muscle.System(this);
            extradofs = 6+56;
            if this.HasSpindle
                extradofs = extradofs + 9;
            end
            extra = extradofs*this.nfibres;
            M = blkdiag(M,speye(extra));
        end
    end
    
end