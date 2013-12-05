function [n,i,t] = getRateParams()
%% intro
% The model uses 2 parameters lists, one per each module
%   This allows for separate simultion of an individual module
%     if you want to do so.
% return: np and ip

 n = zeros(62,1); % NFkB Activation Module
 i = zeros(52,1); % IKK Activation Module
 t = zeros(9,1); %TNF production module 

%% ----- IkB:NFkB Module (including A20 synthesis/deg) -----
%           ikba       ikbb        ikbe
n(1:3)  = [7e-5;      1e-5;       1e-6 ]; % constitutive txn
n(4:6)	= [8;         0.02;       0.3  ]; % inducible txn
n(7:9)  = [3;         3;          3    ]; % Hill Coefficient
n(10:12)= [0;         37;         37   ]; % inducible txn delay
n(13:15)= [.035;      3e-3;       4e-3 ]; % mRNA Degradation
n(16:18)= [0.25;      0.25;       0.25 ]; % translation rate

n(19:21)= [0.09;      9e-3;       0.045]; % Free IkB Import
n(23:25)= [0.012;     0.012;      0.012]; % Free IkB Export
n(27:29)= [0.276;     0.0276;     0.138]; % IkB:NFkB Import
n(30:32)= [0.828;     0.414;      0.414]; % IkB:NFkB Export
n(22)   = 5.4;   % Free NFkB import
n(26)   = 0.0048;% Free NFkB export

n(33:35)= [0.12;      0.18;       0.18 ]; % IkB deg cytoplasm
n(36:38)= [0.12;      0.18;       0.18 ]; % IkB deg nucleus
n(39:41)= [6e-5;      6e-5;       6e-5 ]; % Bound IkB deg cyt
n(42:44)= [6e-5;      6e-5;       6e-5 ]; % Bound IkB deg nuc

n(45:47)= [30;        30;         30   ]; % IkB:NFkB Asn cyt
n(48:50)= [30;        30;         30   ]; % IkB:NFkB Asn nuc
n(51:53)= [6e-5;      6e-5;       6e-5 ]; % IkB:NFkB Dsn cyt
n(54:56)= [6e-5;      6e-5;       6e-5 ]; % IkB:NFkB Dsn nuc

n(57:59)= [0.36;      0.12;       0.18 ]; % Free IkB + IKK deg
n(60:62)= [0.36;      0.12;       0.18 ]; % IkB:NFkB + IKK deg

% A20 (NEW)
n(63)  = 2e-6;    % basal txn rate 2e-6
n(64)  = 0.4;     % inducible txn for A(20),0.4
n(65)  = 3;       % Hill Coefficient
n(66)  = 0;       % inducible txn delay
n(67)  = 0.035;   % mRNA Degradation
n(68)  = 0.25;    % translation rate
n(69)  = 30;      % tsl delay
n(70)  = 0.0029;  % protein degradation
n(71)  = 120;     % promoter shutdown (experiments show ~120min)

%% ----- IKK Activation Module -----
i(1) =     0.1698;  % k_b_lps 1
i(2) =    0.17765;  % k_in_lps 3 
i(3) =    0.26114 ;  % k_out_lps 4
i(4) =     13.367;  % k_d_lpsen, 6

% ligand binding,7:10
i(5)   = 0.19 ;   % k_b_tlr4lps, 7
i(6)   = 2.7  ;   % k_ub_tlr4lps, 8 
i(7)   = 0.19 ;   % k_b_tlr4lpsen, 9
i(8)   = 2.7  ;   % k_ub_tlr4lpsen, 10 

% generation and degradation
i(9)   =   0.025525 ;   % k_g_tlr4, 11
i(10)  =    0.89603;   % k_d_tlr4,12
i(11)  =     2.9295;   % k_d_tlr4en, 13

% shutting
i(12)  =     0.1344  ;   % k_in_tlr4, 14
i(13)  =     3.6099  ;   % k_out_tlf4, 15
i(14)  =   0.23513  ;     % k_in_tlr4lps, 16
i(15)  =   0.041496  ;   % k_out_tlr4lps or 0.02, 17
i(16)  =     14.414;     %k_deg_TLR4LPS 18
i(17)  =    0.42132;    % k_d_tlr4lpsen 19

% activation
i(18)  =   3.2907 ;    % k_f_myd88 TL activated MyD88 -> (MyD88)_6, 20
i(19)  = 3; %Hill coefficient for MyD88 activation 
i(20)  =   0.057936;      %  km_a_myd88, 37
i(21)  =    0.27504;     % k_dis_myd88, 0.02, 21
i(22)  =    0.38509  ;     % k_a_trif  1e+8, 22
i(23)  =   0.011822   ;     % k_i_trif, 0.01, 23

% TRAF6
i(24)  =     4.9488 ; % k_a_TRAF6_MyS88s 24
i(25)  =     1.2988; % k_a_TRAF6_TRIFs 25
i(26)  =    0.21778;       % k_i_TRAF6  27 
                       %i(28)  = 0;         % k_i_TRAF6_A20
% IKKK
i(27)  =    0.34317;   % IKKK_off --> IKKK (TRAF6s mediated), 29
i(28)  =    5e-7;   % IKKK_off --> IKKK (constitutive),30
i(29)  = 0.25;    % IKKK     --> IKKK_off (constitutive), 31

% IKK
i(30)  = 520;    % IKK_off  --> IKK (IKKK mediated),32
i(31)  = 5e-5 ;    % IKK_off  --> IKK (constitutive),33
i(32)  = 0.02  ;    % IKK      --> IKK_off (constitutive), 34
i(33)  = 0.15 ;    % IKK      --> IKK_i (constitutive)  , 35
i(34)  = 0.02  ;    % IKKi     --> IKK_off (constitutive) , 36

%% ----- IRF3  module 
%% IRF3 part
i(35) = 0.91476      ;        %k_a_TBK1_TRIFs, 40 
i(36) = 0.13745    ;      %k_i_TBK1, 41
i(37) = 1.9104    ;      %k_a_IRF_TBK1s, 43
i(38) = 0.00071777        ;      %k_i_IRF,44
i(39) = 0.15531        ;      %k_in_IRFs, 49
i(40) = 2.3393e-07    ;      % k_out_IRFs, 50 
i(41) = 0.002206     ;      %k_i_IRFn, 46
i(42) =      14.85      ;      %k_out_IRF, 48
i(43) =    0.92222      ;      %k_in_IRF, 47
i(44) =    0.40995;     % k_g_irf3, 54
i(45) = 9.0755e-05; % k_d_irf3,
i(46) =  0.0056593; % k_d_irf3s,
i(47) =  0.0028974; % k_d_irf3n,
i(48) =    0.19835; % k_d_irf3ns, 55

%% wt paramters
i(49) =         7.477;      %k_a_TRAF6_MyS88s, 52
i(50) =     3.4132;      %k_a_TRAF6_TRIFs, 53
i(51) =    0.06791; % ikk devider, 38 
i(52) =     633.55; % LPS multiplier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the TNFR parameters
%----- IKK Activation Module ----- NEW f
i(53)   = 0.0154; % pd_m_tnf 45' half life of exogenous TNF

% tnfrm metabolism (synthesis and degradation)
i(54)   = 2e-7;   % tnfrm synthesis (txn, tsl, localization)
i(55)   = 0.0058; % tnfrm --> deg  -- 120' halflife

% TNF-Independent C1 Activation
i(56)   = 1e-5;   % 3tnfrm --> TNFR
i(57)   = 0.1;    % TNFR   --> 3tnfrm
i(58)   = 0.023;  % TNFR internalization -- 30' halflife

i(59)   = 100;    % TNFR + TTR --> C1_off
i(60)   = 0.1;    % C1_off --> TNFR + TTR
i(61)   = 30;     % C1_off --> C1
i(62)  = 2;      % C1     --> C1_off
i(63)  = 1000;   % C1     --> C1_off (A20 Mediated)
i(64)  = i(60);   % C1     --> TNFR + TTR
i(65)  = i(58);   % C1_off internalization
i(66)  = i(58);   % C1 internalization

% TNF-dependent C1 Activation
i(67)  = 1100;   % 3tnfrm + tnf --> TNFRtnf
i(68)  = i(67);  % TNFR + tnf --> TNFRtnf
i(69)  = 0.021;  % TNFRtnf   --> TNFR + tnf
i(70)  = i(58);   % TNFRtnf internalization -- 30' halflife

i(71)  = i(59);   % TNFRtnf + TTR --> C1tnf_off
i(72)  = i(60);   % C1tnf_off --> TNFRtnf + TTR
i(73)  = i(61);   % C1tnf_off --> C1tnf
i(74)  = i(62);  % C1tnf     --> C1tnf_off
i(75)  = i(63);  % C1tnf     --> C1tnf_off (A20 Mediated)
i(76)  = i(60);   % C1tnf     --> TNFRtnf + TTR
i(77)  = i(58);   % C1tnf_off internalization
i(78)  = i(58);   % C1tnf internalization

i(79)  = i(69);  % C1tnf_off --> C1_off + tnf
i(80)  = i(67);  % C1_off + tnf --> C1tnf_off
i(81)  = i(69);  % C1tnf    --> C1 + tnf
i(82)  = i(67);  % C1 + tnf --> C1tnf

% IKKK
i(83)  = 500;    % IKKK_off --> IKKK (C1 mediated)?500
i(84)  = i(83);  % IKKK_off --> IKKK (C1tnf mediated)

%% the CPG-TLR9 parameters
% TLR9rm metabolism (synthesis and degradation)
i(85)   = 2e-7;   % TLR9 synthesis (txn, tsl, localization)
i(86)   = 0.0058; % TLR9 --> deg  -- 120' halflife
i(87)   = .6;     % CpG +  TLR9 -> CpGTLR9, affinity : 186 nM ± 35 nM
i(88)   = 0.1116; % CpGTLR9 -> CpG +  TLR9 1 × 104 M-1 s-1, 6 mins half-life
i(89)   = i(18)*1000;  %MyD88 activation. flux_a_MyD88       = v.IP(20)   * (TLR4LPS)^3 ...
i(90)   = i(20)/100;  %    /((TLR4LPS)^3 + (v.IP(37))^3)  * MyD88;

%% the PIC-TLR3 parameters
% TLR9rm metabolism (synthesis and degradation)
i(91)   = 2e-7;   % TLR3 synthesis (txn, tsl, localization)
i(92)   = 0.0058; % TLR3 --> deg  -- 120' halflife
i(93)   = .6;     % PIC +  TLR3 -> CpGTLR3, affinity : 186 nM ± 35 nM
i(94)   = 0.1116; % CpGTLR3 -> CpG +  TLR3 1 × 104 M-1 s-1, 6 mins half-life
i(95)   = i(22);  % TRIF activation. flux_a_MyD88       = v.IP(20)   * (TLR4LPS)^3 ...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TNF part
%  tnf part (NEW)
t(1) = 0;%1e-5;%tnf constitutive txn 
t(2) = 0;%10;%tnf induced txn 
t(3) = 3;%tnf transcription,Hill coefficient
t(4) = .5;%tnf transcription induction EC50
t(5) = .02;%tnf transcript deg. max rate ~~~ Could be determined by measurement. 
t(6) = .02*5;%tnf nascent process rate 
t(7) = 1;%tnf protein syns rate 
t(8) = .18;%tnf deg rate 
t(9) = .18;%tnf secretion rate 

end
