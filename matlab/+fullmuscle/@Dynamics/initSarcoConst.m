function initSarcoConst(this)
% initSarcoConst: Init method for sacromere constants. 
%
% There are 3 sets of constants, one that is independent from parameter
% `\mu` (representing fibre type) and, for those which depend on `\mu`, one for slow and fast
% twitch fibre parameters respectively.
%
% @author Daniel Wirtz @date 2013-01-09
%
% @new{0,7,dw,2013-01-09} Added this function.
%
% 
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

%% Base constant vector (with fast twitch values from original shorten_org)
c = zeros(105,1);
c(1) = 1.0;
c(2) = 4.8;
c(3) = 150;
c(4) = 0.000001;
c(5) = 0.0025;
c(6) = 0.0005;
c(7) = 96485;
c(8) = 350;
c(9) = 350;
c(10) = 0.0032;
c(11) = 21875;
c(12) = 21875;
c(13) = 1.02;
c(14) = -1.29;
c(15) = 0.0081;
c(16) = 0.288;
c(17) = 0.0131;
c(18) = 4.38;
c(19) = 1.38;
c(20) = 0.067;
c(21) = -46;
c(22) = -40;
c(23) = -45;
c(24) = 70;
c(25) = -78;
c(26) = -40;
c(27) = 150;
c(28) = 5.8;
c(29) = 7.5;
c(30) = 14.7;
c(31) = 9;
c(32) = 10;
c(33) = 7;
c(34) = 18;
c(35) = 40;
c(36) = 8314.41;
c(37) = 293;
c(38) = 19.65;
c(39) = 64.8;
c(40) = 804;
c(41) = 11.1;
c(42) = 0.4;
c(43) = 950;
c(44) = 1;
c(45) = 1;
c(46) = 13;
c(47) = 10;
c(48) = 0.000621;
c(49) = 90;
c(50) = 0.1;
c(51) = 1.0;
c(52) = 0.45;
c(53) = 0.1;
c(54) = 0.1;
c(55) = 0.002;
c(56) = 1000;
c(57) = 0.2;
c(58) = 0.2;
c(59) = 4.5;
c(60) = -20;
c(61) = 4.875;
c(62) = 1;
c(63) = 0.00002;
c(64) = 0.75;
c(65) = 0.75;
c(66) = 1.1;
c(67) = 0.5;
c(68) = 0.04425;
c(69) = 0.115;
c(70) = 140;
c(71) = 0.0417;
c(72) = 0.0005;
c(73) = 1500;
c(74) = 0.000033;
c(75) = 0.003;
c(76) = 0.000004;
c(77) = 0.005;
c(78) = 31000;
c(79) = 0.15;
c(80) = 30;
c(81) = 0.0015;
c(82) = 0.15;
c(83) = 0.375;
c(84) = 1.5;
c(85) = 0;
c(86) = 0.15;
c(87) = 0.15;
c(88) = 0.05;
c(89) = 1.5;
c(90) = 15;
c(91) = 0.24;
c(92) = 0.18;
c(93) = 0.12;
c(94) = 0.00002867;
c(95) = 0.00000362;
c(96) = 1;
c(97) = 0.0001;
c(98) = 6;
c(99) = 300;
c(100) =  0.95*c(66)* pi*(c(67)^2);
c(101) =  0.05*c(66)* pi*(c(67)^2);
c(102) =  0.01*c(100);
c(103) =  0.99*c(100);
c(104) =  0.01*c(101);
c(105) =  0.99*c(101);

positions = [1 2 8 9 10 11 12 13 14 25 28 38 39 40 41 48 49 61 63 68 71 72 74 75 89 90 91 92 93 94 99];
slow = [this.System.C_m_slow 2.79 559 559 0.00174 40229.885 40229.885 0.34 -0.43 ...
    -68 7.1 3.275 10.8 134 1.85 0.0001656 70 2.4375 4e-05 0.0885 0 0 0 0 ...
    .05, .5, 0.008, 0.006, 0.004 ... % Human muscle, for mouse use 0.5 5 0.08 0.06 0.04
    3.94e-06 60];
fast = [this.System.C_m_fast 4.8 350 350 0.0032 21875 21875 1.02 -1.29 -78 5.8 ...
    19.65 64.8 804 11.1 0.000621 90 4.875 2e-05 0.04425 0.0417 0.0005 3.3e-05 0.003 ...
    .15, 1.5, 0.024, 0.018, 0.012 ... % Human muscle, for mouse use 1.5 15 0.24 0.18 0.12
    2.867e-05 300];

% Remainders from former version:
% c3(74) = 0.05;     
% c3(75) = 1.E4;     
% c3(76) = 0; 
% c3(77) = 3.4;  
% c3(78) = 1000;

this.SarcoConst_slow = slow';
this.SarcoConst_fast = fast';
this.SarcoConst_dynpos = positions;
this.SarcoConst_base = c;
end