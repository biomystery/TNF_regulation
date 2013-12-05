function [n,i] = getRateParams()
%% intro
% The model uses 2 parameters lists, one per each module
%   This allows for separate simultion of an individual module
%     if you want to do so.
% return: np and ip

n = zeros(84,1); % NFkB Activation Module
i = zeros(85,1); % IKK Activation Module

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

%  tnf part (NEW)
n(63) = 0;%tnf constitutive txn, 1e-5
n(64) = 0;%tnf induced txn,1 
n(65) = 3;%tnf transcription,Hill coefficient
n(66) = .035;%tnf transcript deg. max rate ~~~ Could be determined by measurement.
n(67) = 1;%tnf transcript deg. TRIF EC50
n(68) = 1/5;%tnf transcript deg. MyD88 EC50
n(69) = .25;%tnf translation rate
n(70) = 1e-2;%tnf basal secretion rate
n(71) = 1;%tnf induced secrection rate
n(72) = 0.1;%tnf induced secrection EC50 for TRIF*
n(73) = 3;%tnf secrection hill cofficient. (trimers)
n(74) = .1;%tnf secrection EC50. (trimers)
n(75) = .18;%tnf protein degradation (inside the cell)

% A20 (NEW)
n(76)  = 2e-6;    % basal txn rate 2e-6
n(77)  = 0.4;     % inducible txn for A(20),0.4
n(78)  = 3;       % Hill Coefficient
n(79)  = 0;       % inducible txn delay
n(80)  = 0.035;   % mRNA Degradation
n(81)  = 0.25;    % translation rate
n(82)  = 30;      % tsl delay
n(83)  = 0.0029;  % protein degradation
n(84)  = 360;     % promoter shutdown (experiments show ~120min)

%% ----- IKK Activation Module -----
% ligand binding
i(7)   = 0.19 ;   % k_b_tlr4lps
i(8)   = 2.7  ;   % k_ub_tlr4lps
i(9)   = 0.19 ;   % k_b_tlr4lpsen
i(10)  = 2.7  ;   % k_ub_tlr4lpsen

% activation
i(20)  = 6 ;    % k_f_myd88 TL activated MyD88 -> (MyD88)_6
i(37)  = 1.5e-3  ;     %  km_a_myd88 
i(21)  = 0.2   ;     % k_dis_myd88, 0.02
i(22)  = 500  ;     % k_a_trif  1e+8
i(23)  = 0.2   ;     % k_i_trif, 0.01
%i(5) belongs to the degradation


% shutting
i(3)   = 0.02  ;   % k_in_lps
i(4)   = 0.02 ;    % k_out_lps
i(14)  = 0.01  ;   % k_in_tlr4
i(15)  = 0.1  ;    % k_out_tlf4
i(16)  = 0.05  ;   % k_in_tlr4lps
i(17)  = 0.005   ; % k_out_tlr4lps or 0.02

% generation and degradation
i(11)  = 2e-6 ;   % k_g_tlr4
i(12)  = 0.001;   % k_d_tlr4
i(13)  = 0.001;   % k_d_tlr4en
i(6)   = log(2)/45;  % k_d_lpsen
i(5)   = log(2)/45;  % k_d_lps
i(19)  = 0.1;     % k_d_tlr4lpsen
i(18) = log(2)/60;   %k_deg_TLR4LPS or 0i

% TRAF6
i(26)  = 5e-8 ;     % k_a_TRAF6; NOT USE; no basal activation for TRAF6
i(24)  = 0.1   ;    % k_a_TRAF6_MyS88s %5000
i(25)  = 0.04 ;     % k_a_TRAF6_TRIFs %50
i(27)  = 0.3;       % k_i_TRAF6 %0.02
i(28)  = 0;         % k_i_TRAF6_A20


% IKKK
i(30)  = 1e-7;%5e-8;   % IKKK_off --> IKKK (constitutive)
i(29)  = 0.05;   % IKKK_off --> IKKK (TRAF6s mediated)
i(31)  = .25;%0.2;    % IKKK     --> IKKK_off (constitutive)

% IKK
i(33)  = 1e-6 ;    % IKK_off  --> IKK (constitutive)
i(32)  = 1000;%0.1  ;    % IKK_off  --> IKK (IKKK mediated) 1000
i(34)  = .1;%0.2  ;    % IKK      --> IKK_off (constitutive) 0.02
i(35)  = .15;%0.02 ;    % IKK      --> IKK_i (constitutive)  0.02
i(36)  = .02;%0.2  ;    % IKKi     --> IKK_off (constitutive) 0.2

%% ----- IRF3  module

%% IRF3 part

i(39) = 1e-8     ;      %k_a_TBK1, NOT USE
i(40) = .04     ;        %k_a_TBK1_TRIFs
i(41) = 0.3/10     ;      %k_i_TBK1
i(42) = 1e-6   ;         %k_a_IRF, NOT USE
i(43) = 0.025    ;      %k_a_IRF_TBK1s
i(44) = 0.125     ;      %k_i_IRF
i(45) = 1e-6    ;      %k_a_IRFn, NOT USE
i(46) = 0.125     ;      %k_i_IRFn
i(47) = 0.05      ;      %k_in_IRF
i(48) = 0.05      ;      %k_out_IRF
i(49) = 0.02        ;      %k_in_IRFs
i(50) = 0.0001   ;      % k_out_IRFs
i(51) = 0        ;      %k_i_TBK1_TRAF6    NOT USE

%% additional parameters for the TL degradation
i(1) = 0.072      ;  % k_b_lps


%% IRF degradation and generation

i(54) = 0.0005;     % k_g_irf3
i(55) = log(2)/360; % k_d_irf3ns

%%1005
i([ 40:41 43:44 46:50 52:55])=[
    0.024888
    0.078294
    0.076191
    0.0056022
    0.0069368
    0.14875
    1.6743
    0.037844
    9.919e-06
    0.002575
    0.08462
    0.028263
    0.042255
    ];


%% not used yet
i([2,5,26,39,42,45,51] ) = [0;0;0;0;0;0;0];

i(48) = 1.5;
i([1 3:4 6:8 29:31 11:25 27 32:38])=[
    0.070092
    0.069044
    4.3852e-05
    0.00024703
    0.19903
    3.1149
    0.044635
    7.6517e-05
    0.19586
    5.0499e-05
    0.011257
    0.0053776
    3.6781e-05
    0.00015059
    0.024577
    0.21779
    0.47681
    3.3114
    79.273
    1.7258
    253.59
    0.037704
    0.001092
    0.024491
    0.19769
    0.18315
    0.0016916
    0.053523
    0.13938
    0.34906
    0.0010089
    0.39071
    ];

% IKKK
i(30)  = 1e-7;%5e-8;   % IKKK_off --> IKKK (constitutive)
i(29)  = 0.05;   % IKKK_off --> IKKK (TRAF6s mediated)
i(31)  = .25;%0.2;    % IKKK     --> IKKK_off (constitutive)

% IKK
i(33)  = 1e-6 ;    % IKK_off  --> IKK (constitutive)
i(32)  = 1000;%0.1  ;    % IKK_off  --> IKK (IKKK mediated) 1000
i(34)  = .1;%0.2  ;    % IKK      --> IKK_off (constitutive) 0.02
i(35)  = .15;%0.02 ;    % IKK      --> IKK_i (constitutive)  0.02
i(36)  = .02;%0.2  ;    % IKKi     --> IKK_off (constitutive) 0.2
i(38)  = 1 ; 

%----- IKK Activation Module -----
i(56)   = 0.0154; % pd_m_tnf 45' half life of exogenous TNF

% tnfrm metabolism (synthesis and degradation)
i(57)   = 2e-7;   % tnfrm synthesis (txn, tsl, localization)
i(58)   = 0.0058; % tnfrm --> deg  -- 120' halflife

% TNF-Independent C1 Activation
i(59)   = 1e-5;   % 3tnfrm --> TNFR
i(60)   = 0.1;    % TNFR   --> 3tnfrm
i(61)   = 0.023;  % TNFR internalization -- 30' halflife

i(62)   = 100;    % TNFR + TTR --> C1_off
i(63)   = 0.1;    % C1_off --> TNFR + TTR
i(64)   = 30;     % C1_off --> C1
i(65)  = 2;      % C1     --> C1_off
i(66)  = 1000;   % C1     --> C1_off (A20 Mediated)
i(67)  = i(63);   % C1     --> TNFR + TTR
i(68)  = i(61);   % C1_off internalization
i(69)  = i(61);   % C1 internalization

% TNF-dependent C1 Activation
i(70)  = 1100;   % 3tnfrm + tnf --> TNFRtnf
i(71)  = i(70);  % TNFR + tnf --> TNFRtnf
i(72)  = 0.021;  % TNFRtnf   --> TNFR + tnf
i(73)  = i(61);   % TNFRtnf internalization -- 30' halflife

i(74)  = i(62);   % TNFRtnf + TTR --> C1tnf_off
i(75)  = i(63);   % C1tnf_off --> TNFRtnf + TTR
i(76)  = i(64);   % C1tnf_off --> C1tnf
i(77)  = i(65);  % C1tnf     --> C1tnf_off
i(78)  = i(66);  % C1tnf     --> C1tnf_off (A20 Mediated)
i(79)  = i(63);   % C1tnf     --> TNFRtnf + TTR
i(80)  = i(61);   % C1tnf_off internalization
i(81)  = i(61);   % C1tnf internalization

i(82)  = i(72);  % C1tnf_off --> C1_off + tnf
i(83)  = i(70);  % C1_off + tnf --> C1tnf_off
i(84)  = i(72);  % C1tnf    --> C1 + tnf
i(85)  = i(70);  % C1 + tnf --> C1tnf

% IKKK
% i(31)  = 1e-7;   % IKKK_off --> IKKK (constitutive)
i(86)  = 500;    % IKKK_off --> IKKK (C1 mediated)?500
i(87)  = i(86);  % IKKK_off --> IKKK (C1tnf mediated)
% i(34)  = 0.25;   % IKKK     --> IKKK_off (constitutive)

%% the CPG-TLR9 parameters
% TLR9rm metabolism (synthesis and degradation)
i(88)   = 2e-7;   % TLR9 synthesis (txn, tsl, localization)
i(89)   = 0.0058; % TLR9 --> deg  -- 120' halflife
i(90)   = .6;     % CpG +  TLR9 -> CpGTLR9, affinity : 186 nM ± 35 nM
i(91)   = 0.1116; % CpGTLR9 -> CpG +  TLR9 1 × 104 M-1 s-1, 6 mins half-life
i(92)   = i(20);  %MyD88 activation. flux_a_MyD88       = v.IP(20)   * (TLR4LPS)^3 ...
i(93)   = i(37);  %    /((TLR4LPS)^3 + (v.IP(37))^3)  * MyD88;

%% the PIC-TLR3 parameters
% TLR9rm metabolism (synthesis and degradation)
i(94)   = 2e-7;   % TLR3 synthesis (txn, tsl, localization)
i(95)   = 0.0058; % TLR3 --> deg  -- 120' halflife
i(96)   = .6;     % PIC +  TLR3 -> CpGTLR3, affinity : 186 nM ± 35 nM
i(97)   = 0.1116; % CpGTLR3 -> CpG +  TLR3 1 × 104 M-1 s-1, 6 mins half-life
i(98)   = i(22);  % TRIF activation. flux_a_MyD88       = v.IP(20)   * (TLR4LPS)^3 ...

end




