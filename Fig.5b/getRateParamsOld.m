function [n,i,t] = getRateParams()
%% intro
% The model uses 2 parameters lists, one per each module
%   This allows for separate simultion of an individual module
%     if you want to do so.
% return: np and ip

 n = zeros(75,1); % NFkB Activation Module
 i = zeros(51,1); % IKK Activation Module
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

%  tnf part (NEW)
t(1) = 1e-5;%tnf constitutive txn 
t(2) = 1;%tnf induced txn 
t(3) = 3;%tnf transcription,Hill coefficient
t(4) = .05;%tnf transcription induction EC50
t(5) = .02;%tnf transcript deg. max rate ~~~ Could be determined by measurement. 
t(6) = .02;%tnf nascent process rate 
t(7) = 1;%tnf protein syns rate 
t(8) = .18;%tnf deg rate 
t(9) = .18;%tnf secretion rate 


%% ----- IKK Activation Module -----
% ligand binding
i(7)   = 0.19 ;   % k_b_tlr4lps
i(8)   = 2.7  ;   % k_ub_tlr4lps
i(9)   = 0.19 ;   % k_b_tlr4lpsen
i(10)   = 2.7  ;   % k_ub_tlr4lpsen

% activation
i(20)   = 6 ;    % k_f_myd88 TL activated MyD88 -> (MyD88)_6
i(37)   = 1.5e-3  ;     %  km_a_myd88 // deleted!!
i(21)   = 0.2   ;     % k_dis_myd88, 0.02
i(22)   = 500  ;     % k_a_trif  1e+8
i(23)  = 0.2   ;     % k_i_trif, 0.01
%i(5) belongs to the degradation


% shutting
i(3)  = 0.02  ;   % k_in_lps
i(4)  = 0.02 ;    % k_out_lps
i(14)  = 0.01  ;   % k_in_tlr4
i(15)  = 0.1  ;   % k_out_tlf4
i(16)  = 0.05  ;     % k_in_tlr4lps
i(17)  = 0.005   ;   % k_out_tlr4lps or 0.02

% generation and degradation
i(11)  = 2e-6 ;   % k_g_tlr4
i(12)  = 0.001;   % k_d_tlr4
i(13)  = 0.001;   % k_d_tlr4en
i(6)  = log(2)/45;  % k_d_lpsen
i(5)  = log(2)/45;  % k_d_lps
i(19)  = 0.1;     % k_d_tlr4lpsen

% TRAF6
i(26)  = 5e-8 ;     % k_a_TRAF6
i(24)  = 0.1   ; % k_a_TRAF6_MyS88s %5000
i(25)  = 0.04 ; % k_a_TRAF6_TRIFs %50
i(27)  = 0.3;       % k_i_TRAF6 %0.02
i(28)  = 0;         % k_i_TRAF6_A20


% IKKK
i(30)  = 5e-8;   % IKKK_off --> IKKK (constitutive)
i(29)  = 0.05;   % IKKK_off --> IKKK (TRAF6s mediated)
i(31)  = 0.2;    % IKKK     --> IKKK_off (constitutive)

% IKK
i(33)  = 1e-6 ;    % IKK_off  --> IKK (constitutive)
i(32)  = 0.1  ;    % IKK_off  --> IKK (IKKK mediated) 1000
i(34)  = 0.2  ;    % IKK      --> IKK_off (constitutive) 0.02
i(35)  = 0.02 ;    % IKK      --> IKK_i (constitutive)  0.02
i(36)  = 0.2  ;    % IKKi     --> IKK_off (constitutive) 0.2

%% ----- IRF3  module 

%% IRF3 part

i(39) = 1e-8     ;      %k_a_TBK1
i(40) = .04     ;        %k_a_TBK1_TRIFs
i(41) = 0.3/10     ;      %k_i_TBK1
i(42) = 1e-6   ;         %k_a_IRF
i(43) = 0.025    ;      %k_a_IRF_TBK1s
i(44) = 0.125     ;      %k_i_IRF
i(45) = 1e-6    ;      %k_a_IRFn
i(46) = 0.125     ;      %k_i_IRFn
i(47) = 0.05      ;      %k_in_IRF
i(48) = 0.05      ;      %k_out_IRF
i(49) = 0.02        ;      %k_in_IRFs
i(50) = 0.0001   ;      % k_out_IRFs
i(51) = 0        ;      %k_i_TBK1_TRAF6    %k_out_IRFs

%% additional parameters for the TL degradation
i(18) = log(2)/60;      %k_deg_TLR4LPS or 0i
i(1) = 0.072      ;  % k_b_lps


%% IRF degradation and generation 

i(54) = 0.0005;     % k_g_irf3
i(55) = log(2)/360; % k_d_irf3ns

%%not use
% i(38) =0.38397; % will used in the scaling 
% 52,53 used for wt paramters
i(52:54) = 0; 
i([1:6 11:27 29:37])=[
% i([1:6 11:27 29:37])=[
   0.094633
   3.1229e-09
      0.11118
    0.0017983
   3.2062e-07
   9.6537e-07
   2.3917e-05
    0.0078275
     0.036891
   1.3416e-06
    0.0018652
      0.06919
    0.0032038
      0.52578
       2.6233
       53.216
      0.72374
        466.2
     0.025462
    0.0020936
    0.0098937
   4.9072e-08
      0.22821
     0.019894
   0.00018972
      0.18459
      0.18201
    0.0010117
     0.044108
      0.10471
      0.37053
   0.00074739

];
i(52:53)=[    0.0052644% wt parameters
      0.03844 %      0.52659
];
n(1:6)=[  0.00092737
   8.1998e-07
   2.0908e-06
       131.66
      0.80807
       0.2745
];
i(38) = 0.55987;

%%1005
i([1 3:4 6 11:25 27 29:38 40:41 43:44 46:50 52:55])=[
     0.070032
      0.06984
   4.9603e-05
    0.0010293
   2.9267e-05
      0.01203
    0.0033402
   3.4723e-05
    0.0004929
     0.042158
       0.2016
      0.47553
       3.3991
       77.808
       1.7669
       251.05
      0.03755
   0.00098305
     0.019807
      0.19735
     0.043331
   0.00010207
      0.17995
      0.27952
    0.0015508
     0.054461
      0.12605
       0.3776
   0.00061645
       1.2389
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
n(1:6)=[   3.6661e-07
   5.5832e-05
   2.8511e-07
       49.942
      0.06832
        1.004

];

n([7           61           13           60           16            4           33           57           30           14])=[
     1.6727
    0.0071837
     0.068929
       1.0737
      0.26987
       41.429
       0.7425
      0.49964
        2.027
    0.0080319
 ];
i([11     7    32    25    29    37    22    24    12    16])=[
       5.3904e-05
      0.17772
      0.17938
     0.024489
     0.038272
    0.0010125
       257.24
    0.0010722
     0.011737
     0.021507
];     
%% not used yet 
i([2,5,26,39,42,45,51] ) = [0;0;0;0;0;0;0]; 
n([7,13,14,16,30,33,57,60,61]) = [3;.035;3e-3;0.25;.828;.12;.36;.36;.12]; 
n(1:6)=[ 7.3172e-07
   1.9235e-05
   1.5454e-09
       16.463
    0.0071388
         0.14
    ];

i(7) = 0.19;
n(1:6)=[7e-05
        1e-05
        1e-06
            8
         0.02
          0.3
];
i([38,48]) = [  0.35318;1.5];
i(7:8)=[      0.19865
        3.081];
% % 110512
i([1 3:4 6:8 20:23 29:31 11:25 27 32:38])=[
        0.070092
     0.069044
   4.3852e-05
   0.00024703
      0.19903
       3.1149
         77.808
       1.7636
       257.25
     0.037659
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



end




