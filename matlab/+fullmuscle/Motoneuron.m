classdef Motoneuron < KerMorObject
% Motoneuron: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-09-24
%
% @new{0,7,dw,2014-09-24} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess=private)
        Constants;
        TypeNoise;
        BaseNoise;
        FibreTypeNoiseFactors;
    end
    
    properties(Access=private)
        % Number of types
        nt;
    end
    
    methods
        
%         function this = Motoneuron
%             this.NoiseGen = models.motoneuron.NoiseGenerator;
%         end
        
        function dy = dydt(this, y, ~)
            c = this.Constants;
            dy = zeros(6,this.nt);
            % dendrites
            dy(1,:) = (-c(1,:).*(y(1,:)-c(11,:))-c(5,:).*(y(1,:)-y(2,:)))./c(7,:);
            % soma
            dy(2,:) = (-c(6,:).*(y(2,:)-c(11,:))-c(5,:).*(y(2,:)-y(1,:))...
                -c(4,:).*y(3,:).^3.*y(4,:).*(y(2,:)-c(9,:))...
                -c(2,:).*y(5,:).^4.*(y(2,:)-c(10,:))...
                -c(3,:).*y(6,:).^2.*(y(2,:)-c(10,:)))./c(8,:);
            % the four gating variables
            dy(3,:) = 0.32*(13-y(2,:))./(exp((13-y(2,:))/5)-1).*(1-y(3,:))...
                -0.28*(y(2,:)-40)./(exp((y(2,:)-40)/5)-1).*(y(3,:));
            dy(4,:) = 0.128*(exp((17-y(2,:))/18)).*(1-y(4,:))-4./(exp((40-y(2,:))/5)+1).*(y(4,:));
            dy(5,:) = 0.032*(15-y(2,:))./(exp((15-y(2,:))/5)-1).*(1-y(5,:))...
                -0.5*(exp((10-y(2,:))/40)).*(y(5,:));
            dy(6,:) = 3.5./(exp((55-y(2,:))/4)+1).*(1-y(6,:))-0.025*(y(6,:));
        end
        
        function J = Jdydt(this, y, ~, typeidx)
            c = this.Constants(:,typeidx);
            J = zeros(6,6);    
            J(1,1) = -(c(1) + c(5))/c(7);
            J(2,1) = c(5)/c(8);
            J(1,2) = c(5)/c(7);
            J(2,2) = -(c(4)*y(4)*y(3)^3 + c(2)*y(5)^4 + c(3)*y(6)^2 + c(5) + c(6))/c(8);
            J(3,2) = (8*(y(3) - 1))/(25*(exp(13/5 - y(2)/5) - 1)) - (7*y(3))/(25*(exp(y(2)/5 - 8) - 1)) + (y(3)*exp(y(2)/5 - 8)*((7*y(2))/25 - 56/5))/(5*(exp(y(2)/5 - 8) - 1)^2) + (exp(13/5 - y(2)/5)*((8*y(2))/25 - 104/25)*(y(3) - 1))/(5*(exp(13/5 - y(2)/5) - 1)^2);
            J(4,2) = (8*exp(17/18 - y(2)/18)*(y(4) - 1))/1125 - (4*y(4)*exp(8 - y(2)/5))/(5*(exp(8 - y(2)/5) + 1)^2);
            J(5,2) = (y(5)*exp(1/4 - y(2)/40))/80 + (4*(y(5) - 1))/(125*(exp(3 - y(2)/5) - 1)) + (exp(3 - y(2)/5)*((4*y(2))/125 - 12/25)*(y(5) - 1))/(5*(exp(3 - y(2)/5) - 1)^2);
            J(6,2) = -(7*exp(55/4 - y(2)/4)*(y(6) - 1))/(8*(exp(55/4 - y(2)/4) + 1)^2);
            J(2,3) = (3*c(4)*y(3)^2*y(4)*(c(9) - y(2)))/c(8);
            J(3,3) =  ((8*y(2))/25 - 104/25)/(exp(13/5 - y(2)/5) - 1) - ((7*y(2))/25 - 56/5)/(exp(y(2)/5 - 8) - 1);
            J(2,4) = (c(4)*y(3)^3*(c(9) - y(2)))/c(8);
            J(4,4) =  - (16*exp(17/18 - y(2)/18))/125 - 4/(exp(8 - y(2)/5) + 1);
            J(2,5) = (4*c(2)*y(5)^3*(c(10) - y(2)))/c(8);
            J(5,5) = ((4*y(2))/125 - 12/25)/(exp(3 - y(2)/5) - 1) - exp(1/4 - y(2)/40)/2;
            J(2,6) = (2*c(3)*y(6)*(c(10) - y(2)))/c(8);
            J(6,6) = - 7/(2*(exp(55/4 - y(2)/4) + 1)) - 1/40;
        end
        
        function setType(this, fibretypes)
            % Set the fibre type of the motoneuron
            %
            % Parameters:
            % fibretypes: Values are assumed to be in [0,1] @type
            % rowvec<double>

            %% Assemble constants
            % Membrane capacitance (NOT to be confused with Cm of the
            % sarcomere model!)
            Cm=1;
            
            Ri=70/1000;
            this.nt = length(fibretypes);
            c = zeros(11,this.nt);
            % cf. Cisi and Kohn 2008, Table 2, page 7
            Rmd = 14.4+6.05-this.coolExp(6.05,14.4,fibretypes);
            Rms=1.15+0.65-this.coolExp(0.65,1.15,fibretypes);
            
            ld=this.coolExp(0.55,1.06,fibretypes);
            ls=this.coolExp(77.5e-6*100,113e-6*100,fibretypes);
            
            rd=this.coolExp(41.5e-6*100,92.5e-6*100,fibretypes)/2;
            rs=this.coolExp(77.5e-6*100,113e-6*100,fibretypes)/2;
            
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
            
            this.Constants = c;
            
            %% Assemble noise signal for each fibre
            ng = models.motoneuron.NoiseGenerator;
            ng.setFibreType(fibretypes(1));
            thenoise = zeros(this.nt,length(ng.indepNoise));
            thenoise(1,:) = ng.indepNoise;
            for k=2:this.nt
                ng.setFibreType(fibretypes(k));
                thenoise(k,:) = ng.indepNoise;
            end
            this.TypeNoise = thenoise;
            this.BaseNoise = ng.baseNoise;
            
            % The noise signal added to the soma is multiplied with a
            % fibretype-dependent factor.
            % These factors are made accessible where needed
            this.FibreTypeNoiseFactors = 1./(pi*this.coolExp(77.5e-4, 0.0113, fibretypes).^2);
        end
    end
    
    methods(Access=private)
        function v = coolExp(~, a, b, mu)
            v = exp(log(100)*mu)*(b-a)/100 + a;
        end
    end
    
end