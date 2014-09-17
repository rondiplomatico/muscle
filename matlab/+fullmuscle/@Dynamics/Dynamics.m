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
        motoconst;
        sarcoconst;
        nfibres;
        
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
        
        forces_scaling;
        forces_scaling_poly = [-28.9060   53.8167  -24.1155   -7.2909    7.3932];
    end

    properties(Transient, Access=private)
        
    end
    
    methods
        function this = Dynamics(sys)
            this = this@muscle.Dynamics(sys);
            this.initSarcoConst;
        end
        
        function configUpdated(this)
            sys = this.System;
            ft = sys.Model.FibreTypes;
            this.nfibres = length(ft);
            
            configUpdated@muscle.Dynamics(this);
            
            % fDim and xDim are from muscle.Dynamics, so add moto+sarco
            this.fDim = this.fDim + (6+56)*this.nfibres;
            this.xDim = this.xDim + (6+56)*this.nfibres;
            this.motoconst = this.getMotoConst(ft);
            this.sarcoconst = this.getSarcoConst(ft);
            
            this.moto_sarco_link_moto_out = sys.off_moto + (2:6:6*this.nfibres);
            this.moto_sarco_link_sarco_in = sys.off_sarco + (1:56:56*this.nfibres);
            
            this.forces_scaling = 1./polyval(this.forces_scaling_poly,ft)';
        end
        
        function prepareSimulation(this, mu)
            prepareSimulation@muscle.Dynamics(this, mu);
            % Create function handles for the link function
            diff = this.MSLink_MaxFactor-this.MSLink_MinFactor;
            funstr = sprintf('%g + exp(-(x-%g)^2/150)*%g',...
                this.MSLink_MinFactor,...
                this.MSLink_MaxFactorSignal,diff);
            this.MSLinkFun = eval(['@(x)' funstr ';']);
            funstr = sprintf('-exp(-(x-%g)^2/150)*%g*(x-%g)/75',...
                this.MSLink_MaxFactorSignal,...
                diff,...
                this.MSLink_MaxFactorSignal);
            this.MSLinkFunDeriv = eval(['@(x)' funstr ';']);
        end
        
        function dy = evaluate(this, y, t)
            sys = this.System;
            
            dy = zeros(this.fDim,1);
            
            %% Mechanics
            uvp_pos = 1:sys.num_uvp_dof;
            uvp = y(uvp_pos);
            % Forces from sarcomeres
            forces = y(sys.sarco_output_idx).*this.forces_scaling;
            dy(uvp_pos) = evaluate@muscle.Dynamics(this, [uvp; forces], t);
            
            %% Motoneurons
            moto_pos = sys.off_moto+(1:sys.num_motoneuron_dof);
            ym = reshape(y(moto_pos),6,[]);
            c = this.motoconst;
            dy_m = zeros(6,this.nfibres);
            % dendrites
            dy_m(1,:) = (-c(1,:).*(ym(1,:)-c(11,:))-c(5,:).*(ym(1,:)-ym(2,:)))./c(7,:);
            % soma
            dy_m(2,:) = (-c(6,:).*(ym(2,:)-c(11,:))-c(5,:).*(ym(2,:)-ym(1,:))...
                -c(4,:).*ym(3,:).^3.*ym(4,:).*(ym(2,:)-c(9,:))...
                -c(2,:).*ym(5,:).^4.*(ym(2,:)-c(10,:))...
                -c(3,:).*ym(6,:).^2.*(ym(2,:)-c(10,:)))./c(8,:);
            % the four gating variables
            dy_m(3,:) = 0.32*(13-ym(2,:))./(exp((13-ym(2,:))/5)-1).*(1-ym(3,:))...
                -0.28*(ym(2,:)-40)./(exp((ym(2,:)-40)/5)-1).*(ym(3,:));
            dy_m(4,:) = 0.128*(exp((17-ym(2,:))/18)).*(1-ym(4,:))-4./(exp((40-ym(2,:))/5)+1).*(ym(4,:));
            dy_m(5,:) = 0.032*(15-ym(2,:))./(exp((15-ym(2,:))/5)-1).*(1-ym(5,:))...
                -0.5*(exp((10-ym(2,:))/40)).*(ym(5,:));
            dy_m(6,:) = 3.5./(exp((55-ym(2,:))/4)+1).*(1-ym(6,:))-0.025*(ym(6,:));
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
        end
        
        function J = getStateJacobian(this, y, t)
            sys = this.System;
            
            %% Mechanics
            uvp_pos = 1:sys.num_uvp_dof;
            uvp = [y(uvp_pos); y(sys.sarco_output_idx).*this.forces_scaling];
            J = getStateJacobian@muscle.Dynamics(this, uvp, t);
            
            %% Motoneuron
            J_m = zeros(6,6);
            for k=1:this.nfibres
                moto_pos = sys.off_moto + 6*(k-1) + (1:6);
                ym = y(moto_pos);
                c = this.motoconst(:,k);
                J_m(1,1) = -(c(1) + c(5))/c(7);
                J_m(2,1) = c(5)/c(8);
                J_m(1,2) = c(5)/c(7);
                J_m(2,2) = -(c(4)*ym(4)*ym(3)^3 + c(2)*ym(5)^4 + c(3)*ym(6)^2 + c(5) + c(6))/c(8);
                J_m(3,2) = (8*(ym(3) - 1))/(25*(exp(13/5 - ym(2)/5) - 1)) - (7*ym(3))/(25*(exp(ym(2)/5 - 8) - 1)) + (ym(3)*exp(ym(2)/5 - 8)*((7*ym(2))/25 - 56/5))/(5*(exp(ym(2)/5 - 8) - 1)^2) + (exp(13/5 - ym(2)/5)*((8*ym(2))/25 - 104/25)*(ym(3) - 1))/(5*(exp(13/5 - ym(2)/5) - 1)^2);
                J_m(4,2) = (8*exp(17/18 - ym(2)/18)*(ym(4) - 1))/1125 - (4*ym(4)*exp(8 - ym(2)/5))/(5*(exp(8 - ym(2)/5) + 1)^2);
                J_m(5,2) = (ym(5)*exp(1/4 - ym(2)/40))/80 + (4*(ym(5) - 1))/(125*(exp(3 - ym(2)/5) - 1)) + (exp(3 - ym(2)/5)*((4*ym(2))/125 - 12/25)*(ym(5) - 1))/(5*(exp(3 - ym(2)/5) - 1)^2);
                J_m(6,2) = -(7*exp(55/4 - ym(2)/4)*(ym(6) - 1))/(8*(exp(55/4 - ym(2)/4) + 1)^2);
                J_m(2,3) = (3*c(4)*ym(3)^2*ym(4)*(c(9) - ym(2)))/c(8);
                J_m(3,3) =  ((8*ym(2))/25 - 104/25)/(exp(13/5 - ym(2)/5) - 1) - ((7*ym(2))/25 - 56/5)/(exp(ym(2)/5 - 8) - 1);
                J_m(2,4) = (c(4)*ym(3)^3*(c(9) - ym(2)))/c(8);
                J_m(4,4) =  - (16*exp(17/18 - ym(2)/18))/125 - 4/(exp(8 - ym(2)/5) + 1);
                J_m(2,5) = (4*c(2)*ym(5)^3*(c(10) - ym(2)))/c(8);
                J_m(5,5) = ((4*ym(2))/125 - 12/25)/(exp(3 - ym(2)/5) - 1) - exp(1/4 - ym(2)/40)/2;
                J_m(2,6) = (2*c(3)*ym(6)*(c(10) - ym(2)))/c(8);
                J_m(6,6) = - 7/(2*(exp(55/4 - ym(2)/4) + 1)) - 1/40;
                J = blkdiag(J,J_m);
            end
            
            %% Sarcomeres
            for k=1:this.nfibres
                sarco_pos = sys.off_sarco + 56*(k-1) + (1:56);
                J = blkdiag(J,this.Jac_Sarco(y(sarco_pos),t,this.sarcoconst(:,k)));
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
            s = load(fullfile(fileparts(which(mc.Name)),'SarcoSparsityPattern'));
            J_sarco = s.JP;
            
            for k=1:this.nfibres
                J = blkdiag(J, J_moto);
            end
            for k=1:this.nfibres
                J = blkdiag(J, J_sarco);
            end
            for k=1:this.nfibres
                % first entry of sarco gets 2nd output of motoneuron
                off_sarco = this.System.off_sarco + (k-1)*56 + 1;
                off_moto = this.System.off_moto + (k-1)*6 + 2;
                J(off_sarco,off_moto) = true;
            end
        end
    end
    
    methods(Access=private)
        function c = getMotoConst(this, fibretypes)
            % getMotoConst: private getter function for motoneuron
            % constants, vectorial implementation
            %
            % mu_fibretype values are assumed to be in [0,1]
            sys = this.System;
            
            % Membrane capacitance (NOT to be confused with Cm of the
            % sarcomere model!)
            Cm=1;
            
            Ri=70/1000;
            c = zeros(11,length(fibretypes));
            % cf. Cisi and Kohn 2008, Table 2, page 7
            Rmd = 14.4+6.05-sys.coolExp(6.05,14.4,fibretypes);
            Rms=1.15+0.65-sys.coolExp(0.65,1.15,fibretypes);
            
            ld=sys.coolExp(0.55,1.06,fibretypes);
            ls=sys.coolExp(77.5e-6*100,113e-6*100,fibretypes);
            
            rd=sys.coolExp(41.5e-6*100,92.5e-6*100,fibretypes)/2;
            rs=sys.coolExp(77.5e-6*100,113e-6*100,fibretypes)/2;
            
            c(1,:) = 2*pi*rd.*ld./Rmd;   % para.Gld
            c(2,:) = 4*2*pi*rs.*ls;      % para.Gkf
            c(3,:) = 16*2*pi*rs.*ls;     % para.Gks
            c(4,:) = 30*2*pi*rs.*ls;     % para.Gna
            c(5,:) = 2./(Ri*ld./(pi*rd.^2)+Ri*ls./(pi*rs.^2));     % para.Gc
            c(6,:) = 2*pi*rs.*ls./Rms;   % para.Gls
            c(7,:) = 2*pi*rd.*ld*Cm;     % para.Cd
            c(8,:) = 2*pi*rs.*ls*Cm;     % para.Cs
            s = ones(size(fibretypes));
            c(9,:) = 120*s;     % para.Vna
            c(10,:) = -10*s;     % para.Vk
            c(11,:) = 0*s;     % para.Vl
        end
        
        function sc = getSarcoConst(this, fibretypes)
            sc = repmat(this.SarcoConst_base,1,this.nfibres);
            sc(this.SarcoConst_dynpos,:) = this.SarcoConst_slow*(1-fibretypes) ...
                + this.SarcoConst_fast * fibretypes;
        end
    end
    
end