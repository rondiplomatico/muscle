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
        num_all_dof;
        off_moto;
        off_sarco;
        off_moto_full;
        off_sarco_full;
        
        input_motoneuron_link_idx;
        sarco_output_idx;
        
        % Membrane capacitance. Different values for different fibre types, 
        % due to different action potential propagation speeds
        % C_m is computed parameter-dependent.
        % These constants are for both slow and fast muscle models and are also used in the
        % first entry of the sarcomere constants computed in
        % models.muscle.FibreDynamics.initSarcoConst @type double
        C_m_slow = 0.58;
        C_m_fast = 1;
        
        noiseGen;
        
        % The upper limit polynomial for maximum mean current dependent on
        % the fibre type.
        %
        % Used at the assembly of B to provide suitable coefficient
        % functions.
        %
        % See also: models.motoneuron.experiments.ParamDomainDetection
        upperlimit_poly;
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
            
            ng = models.motoneuron.NoiseGenerator; 
            % TODO
            ng.setFibreType(0);
            this.noiseGen = ng;
            this.Inputs{1} = @ng.getInput;
        end
        
        function configUpdated(this)
            this.nfibres = length(this.Model.Config.FibreTypes);
            
            configUpdated@muscle.System(this);
            
            %% Add input matrix B
            this.B = this.assembleB;
        end
        
        function prepareSimulation(this, mu, inputidx) 
            % Limit mean current depending on fibre type
            mu(2) = min(polyval(this.upperlimit_poly,mu(1)),mu(2));
            
            prepareSimulation@muscle.System(this, mu, inputidx);
        end
        
        function uvwall = includeDirichletValues(this, t, uvw)
            uvwall_mech = includeDirichletValues@muscle.System(this, t, uvw(1:this.num_uvp_dof,:));
            uvwall = [uvwall_mech; uvw(this.num_uvp_dof+1:end,:)];
        end
        
        function [pm, h_mech] = plot(this, t, y, varargin)
            i = inputParser;
            i.KeepUnmatched = true;
            i.addParamValue('PM',[],@(v)isa(v,'PlotManager'));
            %i.addParamValue('F',[]);
            i.parse(varargin{:});
            r = i.Results;
            
            if ~isempty(r.PM)
                pm = r.PM;
            else
                pm = PlotManager(false,2,2);
                pm.LeaveOpen = true;
            end
            varargin(end+1:end+4) = {'Pool', false, 'PM', pm};
            [~, h_mech] = plot@muscle.System(this, t, y, varargin{:});
        end
        
        function v = coolExp(~, a, b, mu)
            v = exp(log(100)*mu)*(b-a)/100 + a;
        end
    end
    
    methods(Access=protected)
        
        function h = initRefinedPlot(this, t, y, r, pm)
            h = [];
            if length(t) > 1
                h = pm.nextPlot('signal','Motoneuron signal','t [ms]','V_m');
                pos = this.num_uvp_glob + (2:6:6*this.nfibres);
                vals = y(pos,:);
                axis(h,[0 t(end) min(vals(:)) max(vals(:))]);
                hold(h,'on');
                
                h2 = pm.nextPlot('force','Action potential','t [ms]','V_m');
                pos = this.num_uvp_glob + this.num_motoneuron_dof + (1:56:56*this.nfibres);
                vals = y(pos,:);
                axis(h2,[0 t(end) min(vals(:)) max(vals(:))]);
                hold(h2,'on');
                
                h3 = pm.nextPlot('force','Activation','t [ms]','A_s');
                pos = this.num_uvp_glob + this.num_motoneuron_dof + (53:56:56*this.nfibres);
                vals = y(pos,:);
                vals = bsxfun(@times, min(1,t), vals);
                axis(h3,[0 t(end) min(vals(:)) max(vals(:))]);
                hold(h3,'on');
                
                h = [h h2 h3];
            end
        end
        
        function refinedPlot(this, h, t, y, r, ts)
            if ts > 1
                time_part = t(1:ts);
                pos = this.num_uvp_glob + (2:6:6*this.nfibres);
                signal = y(pos,1:ts);
                %walpha = mc.FibreTypeWeights(1,:,1) * signal;
                cla(h(1));
                plot(h(1),time_part,signal);
                %plot(h,times,walpha,'LineWidth',2);
                
                pos = this.num_uvp_glob + this.num_motoneuron_dof + (1:56:56*this.nfibres);
                force = y(pos,1:ts);
                %walpha = mc.FibreTypeWeights(1,:,1) * signal;
                cla(h(2));
                plot(h(2),time_part,force);
                
                pos = this.num_uvp_glob + this.num_motoneuron_dof + (53:56:56*this.nfibres);
                force = y(pos,1:ts);
                force = bsxfun(@times, min(1,time_part), force);
                mc = this.Model.Config;
                walpha = mc.FibreTypeWeights(1,:,1) * force;
%                force = [force; walpha]
                cla(h(3));
                plot(h(3),time_part,force,'r',time_part,walpha,'b');
%                 plotyy(h(3),time_part,force,time_part,walpha);
            end
        end
        
        function updateDofNums(this, mc)
            updateDofNums@muscle.System(this, mc);
            
            this.num_motoneuron_dof = 6*this.nfibres;
            % Motoneurons are beginning after mechanics
            this.off_moto = this.num_uvp_dof; 

            % Sarcomeres are beginning after motoneurons
            this.num_sarco_dof = 56*this.nfibres;
            this.off_sarco = this.off_moto + this.num_motoneuron_dof;
            
            this.num_all_dof = this.off_sarco + this.num_sarco_dof;
            
            % Get the positions where the input signal is mapped to the
            % motoneurons
            this.input_motoneuron_link_idx = this.off_moto + (2:6:6*this.nfibres);
            
            this.sarco_output_idx = this.off_sarco + (53:56:56*this.nfibres);
        end
        
        function x0 = assembleX0(this)
            x0 = zeros(this.num_uvp_dof ...
                + this.num_motoneuron_dof + this.num_sarco_dof,1);
            % Get muscle x0
            x0(1:this.num_uvp_dof) = assembleX0@muscle.System(this);
            
            % Load x0 coefficients for moto/sarco system from file
            mc = metaclass(this);
            s = load(fullfile(fileparts(which(mc.Name)),'x0coeff.mat'));
            x0_motorunit = dscomponents.AffineInitialValue;
            m = size(s.coeff,1);
            for k=1:m
                x0_motorunit.addMatrix(sprintf('polyval([%s],mu(1))',...
                    sprintf('%g ',s.coeff(k,:))),full(sparse(k,1,1,m,1)));
            end
            
            ft = this.Model.Config.FibreTypes;
            for k=1:this.nfibres
                x0ms = x0_motorunit.evaluate(ft(k));
                % add moto
                x0(this.off_moto + 6*(k-1) + (1:6)) = x0ms(1:6);
                % add sarco
                x0(this.off_sarco + 56*(k-1) + (1:56)) = x0ms(7:end);
            end
        end
        
        function B = assembleB(this)
            % The divisor in both coefficient functions is the old para.CS
            % value!!
            i = this.input_motoneuron_link_idx;
            s = 1./(pi*this.coolExp(77.5e-4, 0.0113, this.Model.Config.FibreTypes).^2);
            B = dscomponents.AffLinInputConv;
            % Base noise input mapping
            B.addMatrix('1',sparse(i,ones(size(i)),s,this.num_all_dof,2));
            
            % Independent noise input mapping with µ_4 as mean current factor
            % We need to restrict the maximum mean input current to a
            % reasonable value for each fibre type. luckily, we know them
            % already and hence can provide extra matrices with suitable
            % coefficient functions.
            %
            % See also: 
            maxvals = polyval(this.upperlimit_poly,this.Model.Config.FibreTypes);
            for k=1:this.nfibres
                B.addMatrix(sprintf('min(%g,mu(4))',maxvals(k)),...
                    sparse(i(k),2,s(k),this.num_all_dof,2));
            end
        end
        
        function Daff = assembleDampingMatrix(this)
            Daff_mech = assembleDampingMatrix@muscle.System(this);
            
            Daff = dscomponents.AffLinCoreFun(this);
            extra = (6+56)*this.nfibres;
            D = blkdiag(Daff_mech.getMatrix(1),sparse(extra,extra));
            Daff.addMatrix('mu(1)',D);
        end
        
        function M = assembleMassMatrix(this)
            M = assembleMassMatrix@muscle.System(this);
            extra = (6+56)*this.nfibres;
            M = blkdiag(M,speye(extra));
        end
    end
    
end