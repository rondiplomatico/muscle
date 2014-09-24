classdef Dynamics < muscle.Dynamics;
% Dynamics: 
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
    
    properties%(Access=private)
        sarcoconst;
        
        % constants for sarcomere submodel specific for slow twitch fibre
        %
        % See also: Dynamics.initSarcoConst
        SarcoConst_slow;
        
        % constants for sarcomere submodel specific for fast twitch fibre
        %
        % See also: Dynamics.initSarcoConst
        SarcoConst_fast;
        
        % basis constant set for sarcomere submodel
        %
        % See also: Dynamics.initSarcoConst
        SarcoConst_base;
        
        % Positions of type-dependent constants
        %
        % See also: Dynamics.initSarcoConst
        SarcoConst_dynpos;
    end
    
    %% Properties for motoneuron - sarcomere linking
    properties
        % The `V_s` value of the motoneuron input at which the MSLink_MaxFactor should be attained
        %
        % @type double
        MSLink_MaxFactorSignal = 40;
        
        % The maximal factor with which the `V_s` value of the motoneuron should be amplified
        % when added to the sarcomere equations
        %
        % @type double
        MSLink_MaxFactor = 7;
        
        % The minimal factor at which the `V_s` value of the motoneuron should be amplified
        % when added to the sarcomere equations.
        %
        % Between both limits, an exponential Gaussian weight is applied with a certain radius,
        % so that the amplification factor has its minimum value at zero and maximum value at
        % MSLink_MaxFactorSignal.
        %
        % See also FibreDynamics.getLinkFactor
        %
        % @type double
        MSLink_MinFactor = .3;
        
        MSLinkFun;
        MSLinkFunDeriv;
        
        moto_sarco_link_moto_out;
        moto_sarco_link_sarco_in;
        spindle_moto_link_moto_in;
        
        FrequencyDetector;
        
        % The factors with which the primary and secondary affarent of the
        % spindle is multiplied before considered a "mean input current"
        % (then added to the external signal)
        %
        % @type rowvec<double> @default [1 1]*1e-2
        SpindleAffarentWeights = [1 1]*1e-2;
    end

    properties(Access=private)
        % The upper limit for the mean input current fed to the motoneuron
        % soma. as in this current version this is the sum of spindle
        % feedback and external signal, the max value needs to be available
        % here to limit the sum.
        % The external signal is added at an higher level in the ode (B*u
        % component), but yet the sum of (spindle+ext_sig) <
        % max_moto_signal, which is why the external signal is accessed
        % here, too.
        max_moto_signals;
    end
    
    methods
        function this = Dynamics(sys)
            this = this@muscle.Dynamics(sys);
            this.initSarcoConst;
        end
        
        function configUpdated(this)
            sys = this.System;
            mc = sys.Model.Config;
            ft = mc.FibreTypes;
            this.nfibres = length(ft);
            
            configUpdated@muscle.Dynamics(this);
            
            % fDim and xDim are from muscle.Dynamics, so add moto+sarco
            this.fDim = this.fDim + (6+56+9)*this.nfibres;
            this.xDim = this.xDim + (6+56+9)*this.nfibres;
            this.sarcoconst = this.getSarcoConst(ft);
            
            this.moto_sarco_link_moto_out = sys.off_moto + (2:6:6*this.nfibres);
            this.moto_sarco_link_sarco_in = sys.off_sarco + (1:56:56*this.nfibres);
            this.spindle_moto_link_moto_in = this.moto_sarco_link_moto_out;
            
            this.FrequencyDetector = fullmuscle.FrequencyDetector(this.nfibres);
            this.lambda_dot_pos = mc.SpindlePositions;
            this.lambda_dot = zeros(1,this.nfibres);
            this.Jlambda_dot = zeros(this.nfibres,6*mc.PosFE.Geometry.DofsPerElement);
        end
        
        function prepareSimulation(this, mu)
            prepareSimulation@muscle.Dynamics(this, mu);
            % Create function handles for the link function
            diff = this.MSLink_MaxFactor-this.MSLink_MinFactor;
            funstr = sprintf('%g + exp(-(x-%g).^2/150)*%g',...
                this.MSLink_MinFactor,...
                this.MSLink_MaxFactorSignal,diff);
            this.MSLinkFun = eval(['@(x)' funstr ';']);
            funstr = sprintf('-exp(-(x-%g).^2/150).*(%g*(x-%g))/75',...
                this.MSLink_MaxFactorSignal,...
                diff,...
                this.MSLink_MaxFactorSignal);
            this.MSLinkFunDeriv = eval(['@(x)' funstr ';']);
            
            sys = this.System;
            % Register the ODE callback for the frequency integration
            slv = sys.Model.ODESolver;
            if ~isa(slv,'solvers.MLWrapper')
                error('Only programmed to work with ML-builtin solvers so far!');
            end
            fd = this.FrequencyDetector;
            slv.odeopts = odeset(slv.odeopts,...
                'OutputFcn',@(t,y,flag)fd.processSignal(t,y'),...
                'OutputSel',this.moto_sarco_link_moto_out);
            
            this.max_moto_signals = polyval(sys.upperlimit_poly,sys.Model.Config.FibreTypes);
        end
        
        function dy = evaluate(this, y, t)
            sys = this.System;
            
            dy = zeros(this.fDim,1);
            
            %% Mechanics
            uvp_pos = 1:sys.num_uvp_dof;
            % Use uvp as argument and also pass in s (=sarco forces)
            uvps = [y(uvp_pos); max(0,y(sys.sarco_output_idx)-sys.sarco_mech_signal_offset)];
            dy(uvp_pos) = evaluate@muscle.Dynamics(this, uvps, t);
            
            %% Motoneurons
            mo = sys.Motoneuron;
            moto_pos = sys.off_moto+(1:sys.num_motoneuron_dof);
            dy_m = mo.dydt(reshape(y(moto_pos),6,[]),t);
            dy(moto_pos) = dy_m(:);
            
            %% Sacromeres
            sarco_pos = sys.off_sarco + (1:sys.num_sarco_dof);
            ys = reshape(y(sarco_pos),56,[]);
            dys = this.dydt_sarcomere(ys, t);
            dy(sarco_pos) = dys(:);
            
            %% Link of motoneurons to sarcomeres
            moto_out = y(this.moto_sarco_link_moto_out);
            fac = min(this.MSLink_MaxFactor,this.MSLinkFun(moto_out));
            signal = fac.*moto_out./this.sarcoconst(1,:)';
            % Add signal to corresponding locations
            dy(this.moto_sarco_link_sarco_in) = dy(this.moto_sarco_link_sarco_in) + signal;
            
            %% Spindles
            sp = sys.Spindle;
            spindle_pos = sys.off_spindle + (1:sys.num_spindle_dof);
            yspindle = reshape(y(spindle_pos),9,[]);
            % Get motoneuron frequency
            fd = this.FrequencyDetector;
%             this.lambda_dot
            dys = sp.dydt(yspindle,t,fd.Frequency,this.lambda_dot,0);
            dy(spindle_pos) = dys(:);

            %% Link of spindle to motoneuron
            spindle_sig = this.SpindleAffarentWeights*sp.getAfferents(yspindle);
            ext_sig = this.mu(4); % sys.Inputs{1}(t);
            eff_spindle_sig = min(spindle_sig,this.max_moto_signals - ext_sig);
%             sys.noise(:,round(t)+1)'.*eff_spindle_sig
            noise_sig = mo.TypeNoise(:,round(t)+1)'.*eff_spindle_sig.*mo.FibreTypeNoiseFactors;
%             fprintf('Spindle->Neuron: adding %g at dy(%d)\n',noise_sig,this.spindle_moto_link_moto_in);
            dy(this.spindle_moto_link_moto_in) = ...
                dy(this.spindle_moto_link_moto_in) + noise_sig';
        end
        
        function J = getStateJacobian(this, y, t)
%             J = this.getStateJacobianFD(y,t);
%             return;
            sys = this.System;
            
            %% Mechanics
            uvp_pos = 1:sys.num_uvp_dof;
            uvps = [y(uvp_pos); max(0,y(sys.sarco_output_idx)-sys.sarco_mech_signal_offset)];
            J = getStateJacobian@muscle.Dynamics(this, uvps, t);
            
            %% Motoneuron
            mo = sys.Motoneuron;
            for k=1:this.nfibres
                moto_pos = sys.off_moto + 6*(k-1) + (1:6);
                J = blkdiag(J,mo.Jdydt(y(moto_pos),t,k));
            end
            
            %% Sarcomeres
            for k=1:this.nfibres
                sarco_pos = sys.off_sarco + 56*(k-1) + (1:56);
                J = blkdiag(J,this.Jac_Sarco(y(sarco_pos),t,this.sarcoconst(:,k)));
            end
            
            %% Spindles
            sp = sys.Spindle;
            freq = this.FrequencyDetector.Frequency;
            for k=1:this.nfibres
                spindle_pos = sys.off_spindle + 9*(k-1) + (1:9);
                [Jspin, Jspin_Ldot] = sp.Jdydt(y(spindle_pos), t, freq(k), this.lambda_dot(k), 0);
                J = blkdiag(J,Jspin);
                J(spindle_pos,1:sys.num_u_dof+sys.num_v_dof) = Jspin_Ldot'*this.Jlambda_dot(k,:);
            end
            
            %% Motoneuron to Sarcomere coupling
            moto_out = y(this.moto_sarco_link_moto_out);
            fac = min(this.MSLink_MaxFactor,this.MSLinkFun(moto_out));
            dfac = this.MSLinkFunDeriv(moto_out);
            dsignal_dmotoout = (dfac .* moto_out + fac)./this.sarcoconst(1,:)';
            for k=1:this.nfibres
                J(this.moto_sarco_link_sarco_in(k),this.moto_sarco_link_moto_out(k)) = dsignal_dmotoout(k);
            end
            
            %% Sarcomere to mechanics coupling
            % The JS matrix is generated during the computation of the
            % mechanics jacobian, as the element/gauss loop is computed
            % there anyways. its not 100% clean OOP, but works for now.
            J(sys.num_u_dof + (1:sys.num_v_dof), sys.off_sarco+(1:sys.num_sarco_dof)) = this.JS;
            
            %% Motoneuron to Spindle coupling
            % The frequency-detector is giving piecewise constant signals,
            % hence there is no detectable change in that direction.
            
            %% Spindle to Motoneuron coupling
            for k=1:this.nfibres
                spindle_pos = sys.off_spindle + 9*(k-1) + (1:9);
                moto_pos = sys.off_moto + 6*(k-1) + 2;
                
%                 spindle_sig = this.SpindleAffarentWeights*sp.getAfferents(y(spindle_pos));
%                 ext_sig = this.mu(4); % sys.Inputs{1}(t);
%                 eff_spindle_sig = min(spindle_sig,this.max_moto_signals - ext_sig);
                
                daf = this.SpindleAffarentWeights*sp.getAfferentsJacobian(y(spindle_pos));
                noise = mo.TypeNoise(k,round(t)+1)*mo.FibreTypeNoiseFactors(k);
                
                J(moto_pos,spindle_pos) = noise*daf;
            end
        end
    end
    
    methods(Access=protected)
        function J = computeSparsityPattern(this)
            J = computeSparsityPattern@muscle.Dynamics(this);
            
            % Neuro
            i = [1,1,2,2,2,2,2,2,3,3,4,4,5,5,6,6];
            j = [1,2,1,2,3,4,5,6,2,3,2,4,2,5,2,6];
            J_moto = sparse(i,j,true,6,6);
            
            % Sarco
            mc = metaclass(this);
            s = load(fullfile(fileparts(which(mc.Name)),'JSP_Sarco'));
            J_sarco = s.JP;
            
            % Spindle
            JSp = this.System.Spindle.JSparsityPattern;
            
            for k=1:this.nfibres
                J = blkdiag(J, J_moto);
            end
            for k=1:this.nfibres
                J = blkdiag(J, J_sarco);
            end
            for k=1:this.nfibres
                J = blkdiag(J, JSp);
            end
            
            sys = this.System;
            
            % Moto -> Sarco link
            for k=1:this.nfibres
                % first entry of sarco gets 2nd output of motoneuron
                off_sarco = sys.off_sarco + (k-1)*56 + 1;
                off_moto = sys.off_moto + (k-1)*6 + 2;
                J(off_sarco,off_moto) = true;
            end
            
            % Sarco -> Mechanics link
            for k=1:this.nfibres
                J(sys.num_u_dof + (1:sys.num_v_dof),sys.sarco_output_idx(k)) = true;
            end
        end
    end
    
    methods(Access=private)
        
        
        function sc = getSarcoConst(this, fibretypes)
            sc = repmat(this.SarcoConst_base,1,this.nfibres);
            sc(this.SarcoConst_dynpos,:) = this.SarcoConst_slow*(1-fibretypes) ...
                + this.SarcoConst_fast * fibretypes;
        end
    end
    
end