classdef Spindle < KerMorObject
% Spindle: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-09-18
%
% @new{0,7,dw,2014-09-18} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess=private)
        spindleConst;
        JSparsityPattern;
        y0;
    end
    
    methods
        
        function this = Spindle
            this.initSpindleConst;
           
            % Sparsity pattern
            i = [1 3 4 5 2 4 3 6 4 7 5 8 3 10 4 10 11 5 10 11 3 4 5 11];
            j = [1 1 1 1 2 2 3 3 4 4 5 5 6 6 7 7 7 8 8 8 9 9 9 9];
            this.JSparsityPattern = sparse(i,j,true,11,11);
            
            % Initial values
            y = zeros(11,1);
            y(1) = 0; %0.5*0.735;
            y(2) = 0; %0.5*0.735;
            y(3) = 0.000634066218078;
            y(4) = 0.000634066218078;
            y(5) = 0.000634066218078;
            y(6) = 0.020731324015575;
            y(7) = 0.020731324015575;
            y(8) = 0.020731324015575;
            y(9) = 0.95; %length
            y(10) = 0; % output phase value, 'primary_afferent'
            y(11) = 0; % output phase value, 'secondary_afferent'
            this.y0 = y;
        end
        
        function dy = dydt(this, y, t, moto_sig, Ldot, Lddot)
            % dydt_spindle: Function for change rate of spindle y.

            %% Init
            moto_freq = moto_sig;

            C = 0.42;
            C = C + (Ldot >= 0)*.58; % C = 1 if Ldot > 0
            
            c = this.spindleConst;
            dy = zeros(size(y));
            
            sgntmp1 = sign(Ldot-y(3,:)./c(1));
            sgntmp2 = sign(Ldot-y(4,:)./c(31));
            sgntmpc = sign(Ldot-y(5,:)./c(61));
            
            %% Eval

            % f_dyn_bag1
            dy(1,:) = (moto_freq^c(22)./(moto_freq^c(22)+c(21)^c(22)) - y(1,:))/c(20);
            % f_stat_bag2
            dy(2,:) = (moto_freq.^c(52)/(moto_freq^c(52)+c(51)^c(52)) - y(2,:))/c(20);
            % f_stat_chain
            f_stat_chain = moto_freq^c(82)/(moto_freq^c(82)+c(81)^c(82));

            %% Bag1
            b1  = c(4)+c(5)*y(1,:);
            g1 = c(7)*y(1,:);
            
            dy(3,:) = c(1)/c(3)*(...
                C*b1.*sgntmp1.*abs(Ldot-y(3,:)./c(1)).^c(15)...
                .*(y(9,:)-c(17)-y(6,:)./c(1)-c(16))...
                + c(2)*(y(9,:)-c(17)-y(6,:)/c(1)-c(18))...
                + c(3)*Lddot + g1 - y(6,:));
            dy(6,:) = y(3,:);

            %% Bag2
            b2  = c(34)+c(35)*y(1,:)+c(36)*y(2,:);
            g2 = c(37)*y(1,:)+c(38)*y(2,:);
            dy(4,:) = c(31)/c(33)*(...
                C*b2.*sgntmp2.*abs(Ldot-y(4,:)./c(31)).^c(45)...
                .*(y(9,:) - c(47)-y(7,:)/c(31)-c(46))...
                + c(32)*(y(9,:)-c(47)-y(7,:)/c(31)-c(48))...
                + c(33)*Lddot + g2 - y(7,:));
            dy(7,:) = y(4,:);

            %% Chain
            bc  = c(64)+c(65)*y(1,:) + c(66)*f_stat_chain;
            gc = c(67)*y(1,:) + c(68)*f_stat_chain;
            dy(5,:)  = c(61)/c(63)*(...
                C*bc.*sgntmpc.*abs(Ldot-y(5,:)./c(61)).^c(75)...
                .* (y(9,:)-c(77)-y(8,:)./c(61)-c(76)) +c(62)...
                .* (y(9,:)-c(77)-y(8,:)./c(61)-c(78)) +c(63) .* Lddot + gc - y(8,:));
            dy(8,:)  = y(5,:);

            % Velocity
            dy(9,:) = Ldot;

            %% Frequency = current change of phase
            % This is the signal being put back onto the motoneurons.

            % compute algebraic frequency from above stuff
            Aff_pot_prim_bag1  = c(14)*(y(6,:)/c(1) - (c(12)-c(17)));
            Aff_pot_prim_bag2  = c(44)*(y(7,:)/c(31) - (c(42)-c(47)));
            Aff_pot_prim_chain = c(74)*(y(8,:)/c(61) - (c(72)-c(77)));
            sum_bag2_chain = Aff_pot_prim_bag2 + Aff_pot_prim_chain;

            % Primary_afferent
            maxhelp = Aff_pot_prim_bag1 < sum_bag2_chain;
            dy(10,:) = (1-maxhelp) .* Aff_pot_prim_bag1 + maxhelp.*sum_bag2_chain ...
                + 0.156*(maxhelp .* Aff_pot_prim_bag1 + (1-maxhelp).*sum_bag2_chain);

            Aff_pot_sec_bag2  = c(53)*(c(41)*c(49)/c(47)*(y(7,:)/c(31)-(c(42)-c(47)))...
                + (1-c(41))*c(49)/c(48)*(y(9,:)-y(7,:)/c(31)-c(47)-c(43)));
            Aff_pot_sec_chain = c(83)*(c(71)*c(79)/c(77)*(y(8,:)/c(61)-(c(72)-c(77)))...
                + (1-c(71))*c(79)/c(78)*(y(9,:)-y(8,:)/c(61)-c(77)-c(73)));

            % Secondary_afferent
            dy(11,:) = Aff_pot_sec_bag2 + Aff_pot_sec_chain;

            % The spindle frequency is in Hertz (so per second), but the system time span is millisecond.
            dy(10:11,:) = dy(10:11,:)/1000;

        end
        
        function J = Jdydt(this, y, t, moto_freq, Ldot, Lddot)
            % computes the Jacobian for dydt
            %
          
            moto_freq = 3;

            C = 0.42;
            C = C + (Ldot >= 0)*.58; % C = 1 if Ldot > 0
            
            c = this.spindleConst;
            
            sgntmp1 = sign(Ldot-y(3,:)./c(1));
            sgntmp2 = sign(Ldot-y(4,:)./c(31));
            sgntmpc = sign(Ldot-y(5,:)./c(61));
            
            maxhelp = c(14)*(y(6,:)/c(1) - (c(12)-c(17))) < c(44)*(y(7,:)/c(31) - (c(42)-c(47))) + c(74)*(y(8,:)/c(61) - (c(72)-c(77)));
            
            J = sparse(11,11);
            J(1,1) = -1/c(20);
            J(2,2) = -1/c(20);
            J(3,1) = (c(1)*(c(7) - C*c(5)*sgntmp1*abs(Ldot - y(3)/c(1))^c(15)*(c(16) + c(17) - y(9) + y(6)/c(1))))/c(3);
            J(3,3) = (C*c(15)*sgntmp1*abs(Ldot - y(3)/c(1))^(c(15) - 1)*sign(Ldot - y(3)/c(1))*(c(4) + c(5)*y(1))*(c(16) + c(17) - y(9) + y(6)/c(1)))/c(3);
            J(3,6) = -(c(1)*(c(2)/c(1) + (C*sgntmp1*abs(Ldot - y(3)/c(1))^c(15)*(c(4) + c(5)*y(1)))/c(1) + 1))/c(3);
            J(3,9) = (c(1)*(c(2) + C*sgntmp1*abs(Ldot - y(3)/c(1))^c(15)*(c(4) + c(5)*y(1))))/c(3);
            J(4,1) = (c(31)*(c(37) - C*c(35)*sgntmp2*abs(Ldot - y(4)/c(31))^c(45)*(c(46) + c(47) - y(9) + y(7)/c(31))))/c(33);
            J(4,2) = (c(31)*(c(38) - C*c(36)*sgntmp2*abs(Ldot - y(4)/c(31))^c(45)*(c(46) + c(47) - y(9) + y(7)/c(31))))/c(33);
            J(4,4) = (C*c(45)*sgntmp2*abs(Ldot - y(4)/c(31))^(c(45) - 1)*sign(Ldot - y(4)/c(31))*(c(34) + c(35)*y(1) + c(36)*y(2))*(c(46) + c(47) - y(9) + y(7)/c(31)))/c(33);
            J(4,7) = -(c(31)*(c(32)/c(31) + (C*sgntmp2*abs(Ldot - y(4)/c(31))^c(45)*(c(34) + c(35)*y(1) + c(36)*y(2)))/c(31) + 1))/c(33);
            J(4,9) = (c(31)*(c(32) + C*sgntmp2*abs(Ldot - y(4)/c(31))^c(45)*(c(34) + c(35)*y(1) + c(36)*y(2))))/c(33);
            J(5,1) = (c(61)*(c(67) - C*c(65)*sgntmpc*abs(Ldot - y(5)/c(61))^c(75)*(c(76) + c(77) - y(9) + y(8)/c(61))))/c(63);
            J(5,5) = (C*c(75)*sgntmpc*abs(Ldot - y(5)/c(61))^(c(75) - 1)*sign(Ldot - y(5)/c(61))*(c(64) + c(65)*y(1) + (c(66)*moto_freq^c(82))/(c(81)^c(82) + moto_freq^c(82)))*(c(76) + c(77) - y(9) + y(8)/c(61)))/c(63);
            J(5,8) = -(c(61)*(c(62)/c(61) + (C*sgntmpc*abs(Ldot - y(5)/c(61))^c(75)*(c(64) + c(65)*y(1) + (c(66)*moto_freq^c(82))/(c(81)^c(82) + moto_freq^c(82))))/c(61) + 1))/c(63);
            J(5,9) = (c(61)*(c(62) + C*sgntmpc*abs(Ldot - y(5)/c(61))^c(75)*(c(64) + c(65)*y(1) + (c(66)*moto_freq^c(82))/(c(81)^c(82) + moto_freq^c(82)))))/c(63);
            J(6,3) = 1;
            J(7,4) = 1;
            J(8,5) = 1;
            J(10,6) = (39*c(14)*maxhelp)/(250000*c(1)) - (c(14)*(maxhelp - 1))/(1000*c(1));
            J(10,7) = (c(44)*maxhelp)/(1000*c(31)) - (39*c(44)*(maxhelp - 1))/(250000*c(31));
            J(10,8) = (c(74)*maxhelp)/(1000*c(61)) - (39*c(74)*(maxhelp - 1))/(250000*c(61));
            J(11,7) = (c(53)*((c(49)*(c(41) - 1))/(c(31)*c(48)) + (c(41)*c(49))/(c(31)*c(47))))/1000;
            J(11,8) = (c(83)*((c(79)*(c(71) - 1))/(c(61)*c(78)) + (c(71)*c(79))/(c(61)*c(77))))/1000;
            J(11,9) = - (c(49)*c(53)*(c(41) - 1))/(1000*c(48)) - (c(79)*c(83)*(c(71) - 1))/(1000*c(78));
        end
    end
    
    methods(Access=private)
        function initSpindleConst(this)
            % initSpindleConst: Init function for spindle constants
            
            % copied from spindle_whole_2012_10_11.m, line 217 - 291
            c = zeros(1,83);

            % Non-existent for Bag_1: c(6), c(8), c(11), c(13), c(19)
            % c(9) replaced with C as an input parameter
            c(1:11) = [10.4649, 0.15, 0.0002, 0.0605, 0.2592, 0, 0.0289, 0, 0, 0, 0];
            c(12:22) = [0.0423, 0, 20000, 0.3, 0.46, 0.04, 0.76, 0, 0.149, 60, 2];

            % Non-existent for Bag_2: c(35), c(37); c(39) replaced with C as an input parameter
            c(31:42) = [10.4649, 0.15, 0.0002, 0.0822, 0, -0.046, 0, 0.0636, 0, 0, 0.7, 0.0423];
            c(43:53) = [0.89, 10000, 0.3, 0.46, 0.04, 0.76, 0.04, 0.205, 60, 2, 7250];

            % Non-existent for Chain: c(65), c(67),c(80); c(69) replaced with C as an input parameter
            c(61:72) = [10.4649, 0.15, 0.0002, 0.0822, 0, -0.069, 0, 0.0954, 0, 0, 0.7, 0.0423];
            c(73:83) = [0.89, 10000, 0.3, 0.46, 0.04, 0.76, 0.04, 0, 90, 2, 7250];

            this.spindleConst = c;
        end
    end
    
    methods(Static)
        function test_Spindle
            s = fullmuscle.Spindle;
            t = 0:.1:40;
            moto_sig = 1;
            Ldot = 1;
            Lddot = .3;
            opt = odeset('Jacobian',@(y,t)s.Jdydt(t,y,moto_sig,Ldot,Lddot));
            [t,y] = ode15s(@(t,y)s.dydt(y,t,moto_sig,Ldot,Lddot),t,s.y0,opt);
            plot(t,y);
        end
        
        function res = test_Spindle_Jac
            s = fullmuscle.Spindle;
            d = 11;
            n = 100;
            moto_sig = rand;
            Ldot = rand;
            Lddot = rand;
            testx = rand(d,n)*25;
            res = true;
            for k = 1:n
                x = testx(:,k);
                X = repmat(x,1,d);
                dx = ones(d,1)*sqrt(eps(class(x))).*max(abs(x),1);
                DX = sparse(1:d,1:d,dx,d,d);

                % Evaluate makes use of built-in multi-argument evaluation
                Jc = (s.dydt(X+DX,0,moto_sig,Ldot,Lddot)...
                    - repmat(s.dydt(x,0,moto_sig,Ldot,Lddot),1,d))*diag(1./dx);
                J = s.Jdydt(x,0,moto_sig,Ldot,Lddot);    
                diff = abs((Jc-J)./J);
                diff(isnan(diff)) = 0;
                res = res & max(diff(s.JSparsityPattern)) < .005;
            end
        end
        
        function createSpindleJac
            c = sym('c',[83 1]);
            y = sym('y',[11 1]); 
            dy = sym('dy',[11 1]);

            moto_freq = sym('moto_freq','positive');
            C = sym('C','positive');
            Ldot = sym('Ldot');
            Lddot = sym('Lddot');
            % beta = sym('beta');
            % gamma = sym('gamma');
            % f_stat_chain = sym('f_stat_chain');
            % Aff_pot_prim_bag1 = sym('Aff_pot_prim_bag1');
            % Aff_pot_prim_bag2 = sym('Aff_pot_prim_bag2');
            % Aff_pot_prim_chain = sym('Aff_pot_prim_chain');
            % sum_bag2_chain = sym('sum_bag2_chain');
            % Aff_pot_sec_bag2 = sym('Aff_pot_sec_bag2');
            % Aff_pot_sec_chain = sym('Aff_pot_sec_chain');

            % Helper variables to avoid piecewise-stuff in derivatives
            maxhelp = sym('maxhelp','positive');
            sgntmp1 = sym('sgntmp1','positive');
            sgntmp2 = sym('sgntmp2','positive');
            sgntmpc = sym('sgntmpc','positive');

            % Placeholders for conditional evaluations
            % h73 = sym('h73','positive'); % alg(72) > 0;
            % h74 = sym('h74','positive'); % alg(72) <= 0;
            % h75 = sym('h75','positive'); % alg(31) > 0;

            %% Definitions (copy code here as of ""Eval"" section in dydt)
            % f_dyn_bag1
            dy(1,:) = (moto_freq^c(22)./(moto_freq^c(22)+c(21)^c(22)) - y(1,:))/c(20);
            % f_stat_bag2
            dy(2,:) = (moto_freq.^c(52)/(moto_freq^c(52)+c(51)^c(52)) - y(2,:))/c(20);
            % f_stat_chain
            f_stat_chain = moto_freq^c(82)/(moto_freq^c(82)+c(81)^c(82));

            %% Bag1
            b1  = c(4)+c(5)*y(1,:);
            g1 = c(7)*y(1,:);

            dy(3,:) = c(1)/c(3)*(...
                C*b1.*sgntmp1.*abs(Ldot-y(3,:)./c(1))^c(15)...
                .*(y(9,:)-c(17)-y(6,:)./c(1)-c(16))...
                + c(2)*(y(9,:)-c(17)-y(6,:)/c(1)-c(18))...
                + c(3)*Lddot + g1 - y(6,:));
            dy(6,:) = y(3,:);

            %% Bag2
            b2  = c(34)+c(35)*y(1,:)+c(36)*y(2,:);
            g2 = c(37)*y(1,:)+c(38)*y(2,:);
            dy(4,:) = c(31)/c(33)*(...
                C*b2.*sgntmp2.*abs(Ldot-y(4,:)./c(31))^c(45)...
                .*(y(9,:) - c(47)-y(7,:)/c(31)-c(46))...
                + c(32)*(y(9,:)-c(47)-y(7,:)/c(31)-c(48))...
                + c(33)*Lddot + g2 - y(7,:));
            dy(7,:) = y(4,:);

            %% Chain
            bc  = c(64)+c(65)*y(1,:) + c(66)*f_stat_chain;
            gc = c(67)*y(1,:) + c(68)*f_stat_chain;
            dy(5,:)  = c(61)/c(63)*(...
                C*bc.*sgntmpc.*abs(Ldot-y(5,:)./c(61))^c(75)...
                .* (y(9,:)-c(77)-y(8,:)./c(61)-c(76)) +c(62)...
                .* (y(9,:)-c(77)-y(8,:)./c(61)-c(78)) +c(63) .* Lddot + gc - y(8,:));
            dy(8,:)  = y(5,:);

            % Velocity
            dy(9,:) = Ldot;

            %% Frequency = current change of phase
            % This is the signal being put back onto the motoneurons.

            % compute algebraic frequency from above stuff
            Aff_pot_prim_bag1  = c(14)*(y(6,:)/c(1) - (c(12)-c(17)));
            Aff_pot_prim_bag2  = c(44)*(y(7,:)/c(31) - (c(42)-c(47)));
            Aff_pot_prim_chain = c(74)*(y(8,:)/c(61) - (c(72)-c(77)));
            sum_bag2_chain = Aff_pot_prim_bag2 + Aff_pot_prim_chain;

            % Primary_afferent
            % maxhelp = Aff_pot_prim_bag1 < sum_bag2_chain;
            dy(10,:) = (1-maxhelp) * Aff_pot_prim_bag1 + maxhelp*sum_bag2_chain ...
                + 0.156*(maxhelp * Aff_pot_prim_bag1 + (1-maxhelp)*sum_bag2_chain);

            Aff_pot_sec_bag2  = c(53)*(c(41)*c(49)/c(47)*(y(7,:)/c(31)-(c(42)-c(47)))...
                + (1-c(41))*c(49)/c(48)*(y(9,:)-y(7,:)/c(31)-c(47)-c(43)));
            Aff_pot_sec_chain = c(83)*(c(71)*c(79)/c(77)*(y(8,:)/c(61)-(c(72)-c(77)))...
                + (1-c(71))*c(79)/c(78)*(y(9,:)-y(8,:)/c(61)-c(77)-c(73)));

            % Secondary_afferent
            dy(11,:) = Aff_pot_sec_bag2 + Aff_pot_sec_chain;

            % The spindle frequency is in Hertz (so per second), but the system time span is millisecond.
            dy(10:11,:) = dy(10:11,:)/1000;

            %% Create partial derivatives
            JP = sparse(false(11,11));
            f = 1;
            fprintf(f,'J = sparse(11,11);\n');
            for i = 1:11
                for j = 1:11
                    pd = diff(dy(i),y(j));
                    if pd ~= 0
                        JP(i,j) = true;
                        % Convert y14 to y(14) etc, also for c
                        body = regexprep(char(pd), '(c|y)(\d+)', '$1($2)');
                        fprintf(f,'J(%d,%d) = %s;\n',i,j,body);
                    end
                end
            end
            [i,j] = find(JP);
            sprintf('%d ',i)
            sprintf('%d ',j)
        end
    end
    
end