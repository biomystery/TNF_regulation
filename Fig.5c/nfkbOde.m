% check NEW tag

%% 
function delta = nfkbOde(t,x,ode_options,v)
% We ran the ode15s solver with the default options (an empty [])
%  See the ode15s help file for a list of options you can use
    
% Initialize the iteration
    
% Set the persistant variables that handle the delay reactions
    persistent DELAY_NFKB;
    persistent DELAY_A20;  %NEW
    
    if isempty(t)
        % The sim function calls the ode function once without arguments to
        %   reset the persistent variables. Once done, it returns.
        DELAY_NFKB        = cell(3,1);       % Time, concentration, index
        DELAY_NFKB([1 2]) = {zeros(1000,1)}; % increase if getting errors
        DELAY_NFKB(3)     = {1}; % index (starts at 1)
        
        % NEW
        DELAY_A20         = cell(3,1);       % Time, concentration, index
        DELAY_A20([1 2])  = {zeros(1000,1)}; % increase if getting errors
        DELAY_A20(3)      = {1}; % index (starts at 1)

        return;
    end
    
    % Load Previous Concentrations
    [a,b] = size(x);
    delta = zeros(a,b);
    
    %nf-kb part
    IkBa           = x(1);
    IkBan          = x(2);
    IkBaNFkB       = x(3);
    IkBaNFkBn      = x(4);
    IkBat          = x(5);
    IkBb           = x(6);
    IkBbn          = x(7);
    IkBbNFkB       = x(8);
    IkBbNFkBn      = x(9);
    IkBbt          = x(10);
    IkBe           = x(11);
    IkBen          = x(12);
    IkBeNFkB       = x(13);
    IkBeNFkBn      = x(14);
    IkBet          = x(15);
    NFkB           = x(16);
    NFkBn          = x(17);

    
    % TLR4&IKK part
    LPSo           = x(18); % input LPS,41         
    LPS            = x(19); %18
    LPSen          = x(20); %25
    TLR4           = x(21); %19
    TLR4en         = x(22);        %34
    TLR4LPS        = x(23);%20
    TLR4LPSen      = x(24);        %26
    MyD88          = x(25);%21
    MyD88s         = x(26);%22
    TRIF           = x(27); %MyD88_off=0.1-MyD88,23
    TRIFs          = x(28); %so is trif, 24
    TRAF6          = x(29); %32
    TRAF6s         = x(30);%33
    IKKK_off       = x(31); %kept, 28        
    IKKK           = x(32); %kept, 27
    IKK_off        = x(33); %kept, 30        
    IKK            = x(34); %kept, 29
    IKK_i          = x(35); %kept, 31

    %IRF3 part
    TBK1           = x(36); %35
    TBK1s          = x(37); %36
    IRF3           = x(38); %37
    IRF3s          = x(39); %39        
    IRF3n          = x(40);%38
    IRF3ns         = x(41); %40

    %TNF part (NEW)        
    TNFnas        = x(42);
    TNFmRNA       = x(43);
    TNFpro        = x(44);
    TNF           = x(45);
    
    % TNFR part (NEW, total 8)
    tnfrm          = x(46);
    TNFR           = x(47);
    TNFRtnf        = x(48);
    C1             = x(49);
    C1_off         = x(50);
    C1tnf          = x(51);
    C1tnf_off      = x(52);
    TTR            = x(53);
    
    % A20 (NEW) 
    a20            = x(54);
    a20t           = x(55);

    % CpG-TLR9 (NEW)
    CpG            = x(56); 
    TLR9           = x(57);
    CpGTLR9        = x(58); 

    % PIC-TLR3 (NEW)
    PIC            = x(59); 
    TLR3           = x(60);
    PICTLR3        = x(61); 

    %% Set iteration changes (delta) to zero
    delta_IkBa        = 0;
    delta_IkBan       = 0;
    delta_IkBaNFkB    = 0;
    delta_IkBaNFkBn   = 0;
    delta_IkBat       = 0;
    delta_IkBb        = 0;
    delta_IkBbn       = 0;
    delta_IkBbNFkB    = 0;
    delta_IkBbNFkBn   = 0;
    delta_IkBbt       = 0;
    delta_IkBe        = 0;
    delta_IkBen       = 0;
    delta_IkBeNFkB    = 0;
    delta_IkBeNFkBn   = 0;
    delta_IkBet       = 0;
    delta_NFkB        = 0;
    delta_NFkBn       = 0;
    %TLR4 module
    delta_LPS         = 0;
    delta_TLR4        = 0;
    delta_TLR4LPS     = 0;
    delta_MyD88       = 0;
    delta_MyD88s      = 0;
    delta_TRIF        = 0;
    delta_TRIFs       = 0;
    delta_LPSen       = 0;
    delta_TLR4LPSen   = 0;
    %IKK module
    delta_IKKK        = 0;
    delta_IKKK_off    = 0;
    delta_IKK         = 0;
    delta_IKK_off     = 0;
    delta_IKK_i       = 0;
    %TLR4 continue
    delta_TRAF6       = 0;
    delta_TRAF6s      = 0;
    delta_TLR4en      = 0;
    %IRF module
    delta_TBK1        = 0;
    delta_TBK1s          = 0;
    delta_IRF3           = 0;
    delta_IRF3n          = 0;
    delta_IRF3s          = 0;
    delta_IRF3ns         = 0;
    % LPS binding to the surface
    delta_LPSo           = 0; 

    % TNF component (NEW)
    delta_TNFnas  = 0 ;
    delta_TNFmRNA = 0 ;
    delta_TNFpro  = 0 ;
    delta_TNF  = 0 ;        
    
    % TNFR part (NEW)
    delta_tnfrm       = 0;
    delta_TNFR        = 0;
    delta_TNFRtnf     = 0;
    delta_C1          = 0;
    delta_C1_off      = 0;
    delta_C1tnf       = 0;
    delta_C1tnf_off   = 0;
    delta_TTR         = 0;

    % A20 NEW
    delta_a20         = 0;
    delta_a20t        = 0;

    % CpG-TLR9 (NEW)
    delta_CpG         = 0;
    delta_TLR9        = 0;
    delta_CpGTLR9     = 0;

    % PIC-TLR3 (NEW)
    delta_PIC         = 0;
    delta_TLR3        = 0;
    delta_PICTLR3     = 0;


    
    %% init delay
    % Calculate Inducible Transcription and A20 Translation Delays
    % By default, set the delays to the current values (no delay)
    delayed_nfkbn_a   = NFkBn;
    delayed_nfkbn_b   = NFkBn;
    delayed_nfkbn_e   = NFkBn;
    delayed_nfkbn_a20 = NFkBn; %NEW
    delayed_a20t      = a20t ; %NEW

    if v.PHASE == 2 % no delay in phase 1
        
        % Create a cache of previous NFkBn values for txn delays.
        %   The find command is needed to prevent duplicate values that
        %   break the interpolation function
        if (isempty(find(DELAY_NFKB{1}(1:DELAY_NFKB{3})== t,1)))
            DELAY_NFKB{1}(DELAY_NFKB{3}) = t;
            DELAY_NFKB{2}(DELAY_NFKB{3}) = NFkBn;
            DELAY_NFKB{3}= DELAY_NFKB{3} + 1; % iterate the index + 1
        end
        
        if v.NP(10) > 0  % IkBa inducible txn delay
            if t > v.NP(10)
                delayed_nfkbn_a = ...
                    interp1(DELAY_NFKB{1}(1:DELAY_NFKB{3}-1),...
                            DELAY_NFKB{2}(1:DELAY_NFKB{3}-1),t-v.NP(10));
                % because iterate inde
            else
                delayed_nfkbn_a = ...
                    DELAY_NFKB{2}(1); % 1st index is the basal state value
            end
        end
        
        if v.NP(11) > 0  % IkBb inducible txn delay
            if t > v.NP(11)
                delayed_nfkbn_b = ...
                    interp1(DELAY_NFKB{1}(1:DELAY_NFKB{3}-1),...
                            DELAY_NFKB{2}(1:DELAY_NFKB{3}-1),t-v.NP(11));
            else
                delayed_nfkbn_b = ...
                    DELAY_NFKB{2}(1); % 1st index is the basal state value
            end
        end
        
        if v.NP(12) > 0  % IkBe inducible txn delay
            if t > v.NP(12)
                if v.NP(12) == v.NP(11)  % saves cpu time if IkBb=IkBe delay
                    delayed_nfkbn_e = delayed_nfkbn_b;
                else
                    delayed_nfkbn_e = ...
                        interp1(DELAY_NFKB{1}(1:DELAY_NFKB{3}-1),...
                                DELAY_NFKB{2}(1:DELAY_NFKB{3}-1),t-v.NP(12));
                end
            else
                delayed_nfkbn_e = ...
                    DELAY_NFKB{2}(1); % 1st index is the basal state value
            end
        end
    end
    
    
    if v.NP(66) > 0  % A20 inducible txn delay NEW
        if t > v.NP(66)
            delayed_nfkbn_a20 = ...
                interp1(DELAY_NFKB{1}(1:DELAY_NFKB{3}-1),...
                        DELAY_NFKB{2}(1:DELAY_NFKB{3}-1),t-v.NP(66));
        else
            delayed_nfkbn_a20 = ...
                DELAY_NFKB{2}(1); % 1st index is the basal state value
        end
    end
    
    % This code enables the A20 translational delay (NEW)
    %   It does so by caching A20 mRNA levels (NEW)
    if v.NP([64,68,69]) > 0
        % Create a cache of previous A20 mRNA values.
        if(isempty(find(DELAY_A20{1}(1:DELAY_A20{3}) == t,1)))
            DELAY_A20{1}(DELAY_A20{3}) = t;
            DELAY_A20{2}(DELAY_A20{3}) = delayed_a20t;
            DELAY_A20{3}= DELAY_A20{3} + 1; % iterate the index + 1
        end
        
        % Calculate the A20 mRNA value at time = t-delay
        if t > v.NP(69)
            delayed_a20t = ...
                interp1(DELAY_A20{1}(1:DELAY_A20{3}-1),...
                        DELAY_A20{2}(1:DELAY_A20{3}-1),t-v.NP(69));
        else
            delayed_a20t = DELAY_A20{2}(1);
        end
    end
    

    %% set IKK flux
    if(t>0 && strcmp(v.STIMULI,'LPS'))
        IKK_flux = IKK/(v.IKK_TOTAL*v.IP(51)); % use a unused parameters
    else 
        IKK_flux = IKK/v.IKK_TOTAL; 
    end

    % NEW for the TNF production module 
    if(t>=0 && t<=120)
        kdeg_real = interp1([0 30 60 120],[0.07 v.TP(5) 0.07 ...
                            0.07],t);
    else
        kdeg_real = .07; 
    end


    %% Calculate Iteration Fluxes

    %-- NFkB Activation Module
    flux_rsu_a         = v.NP(1) ;
    flux_rsu_b         = v.NP(2) ;
    flux_rsu_e         = v.NP(3) ;
    flux_rsu_f         = v.TP(1) ; %tnf basal transcription (NEW)

    flux_rsr_an        = v.NP(4)   * (delayed_nfkbn_a ^ v.NP(7)) ;
    flux_rsr_bn        = v.NP(5)   * (delayed_nfkbn_b ^ v.NP(8)) ;
    flux_rsr_en        = v.NP(6)   * (delayed_nfkbn_e ^ v.NP(9)) ;

    flux_rsr_fn        = v.TP(2)  * (delayed_nfkbn_a ^ v.TP(3))/ ...
        (delayed_nfkbn_a ^ v.TP(3) + v.TP(4)^v.TP(3)); %tnf transcription (NEW
                                                       %Same
                                                       %no
                                                       %delay
    flux_rp_f         = v.TP(6) * TNFnas; % NEW tnf RNA processing


    flux_rd_a          = v.NP(13)  * IkBat	;
    flux_rd_b          = v.NP(14)  * IkBbt	;
    flux_rd_e          = v.NP(15)  * IkBet	;
    flux_rd_f          = kdeg_real * TNFmRNA; % NEW tnf mRNA
                                              % degradation 

    flux_ps_c_a        = v.NP(16)	* IkBat	;
    flux_ps_c_b        = v.NP(17)	* IkBbt	;
    flux_ps_c_e        = v.NP(18)	* IkBet	;
    flux_ps_c_f        = v.TP(7)	* TNFmRNA;  %tnf protein synthsis (NEW)
    flux_pd_f          = v.TP(8)    * TNFpro;  % tnf protein deg. (only inside the cell,NEW)
    flux_sc_c_f        = v.TP(9)    * TNFpro;  %tnf protein secrection (NEW)


    flux_in_a          = v.NP(19)	* IkBa	;
    flux_in_b          = v.NP(20)	* IkBb	;
    flux_in_e          = v.NP(21)	* IkBe	;
    flux_in_n          = v.NP(22)	* NFkB	;
    flux_ex_a          = v.NP(23)	* IkBan	;
    flux_ex_b          = v.NP(24)	* IkBbn	;
    flux_ex_e          = v.NP(25)	* IkBen	;
    flux_ex_n          = v.NP(26)	* NFkBn	;
    flux_in_2an        = v.NP(27)	* IkBaNFkB	;
    flux_in_2bn        = v.NP(28)	* IkBbNFkB	;
    flux_in_2en        = v.NP(29)	* IkBeNFkB	;
    flux_ex_2an        = v.NP(30)	* IkBaNFkBn	;
    flux_ex_2bn        = v.NP(31)	* IkBbNFkBn	;
    flux_ex_2en        = v.NP(32)	* IkBeNFkBn	;

    flux_pd_c_a        = v.NP(33)	* IkBa	;
    flux_pd_c_b        = v.NP(34)	* IkBb	;
    flux_pd_c_e        = v.NP(35)	* IkBe	;
    flux_pd_n_a        = v.NP(36)	* IkBan	;
    flux_pd_n_b        = v.NP(37)	* IkBbn	;
    flux_pd_n_e        = v.NP(38)	* IkBen	;
    flux_pd_c_2an      = v.NP(39)	* IkBaNFkB	;
    flux_pd_c_2bn      = v.NP(40)	* IkBbNFkB	;
    flux_pd_c_2en      = v.NP(41)	* IkBeNFkB	;
    flux_pd_n_2an      = v.NP(42)	* IkBaNFkBn	;
    flux_pd_n_2bn      = v.NP(43)	* IkBbNFkBn	;
    flux_pd_n_2en      = v.NP(44)	* IkBeNFkBn	;

    flux_a_c_an        = v.NP(45)  * IkBa	* NFkB	;
    flux_a_c_bn        = v.NP(46)  * IkBb	* NFkB	;
    flux_a_c_en        = v.NP(47)  * IkBe	* NFkB	;
    flux_a_n_an        = v.NP(48)	* IkBan	* NFkBn	;
    flux_a_n_bn        = v.NP(49)	* IkBbn	* NFkBn	;
    flux_a_n_en        = v.NP(50)  * IkBen	* NFkBn	;
    flux_d_c_an        = v.NP(51)	* IkBaNFkB	;
    flux_d_c_bn        = v.NP(52)	* IkBbNFkB	;
    flux_d_c_en        = v.NP(53)	* IkBeNFkB	;
    flux_d_n_an        = v.NP(54)	* IkBaNFkBn	;
    flux_d_n_bn        = v.NP(55)	* IkBbNFkBn	;
    flux_d_n_en        = v.NP(56)	* IkBeNFkBn	;

    % IKK Mediated IkB Degradation (free and bound)
    flux_ph_c_a        = v.NP(57)  * IkBa * IKK_flux ;
    flux_ph_c_b        = v.NP(58)	* IkBb * IKK_flux ;
    flux_ph_c_e        = v.NP(59)	* IkBe * IKK_flux ;

    flux_ph_c_an       = v.NP(60)  * IkBaNFkB * IKK_flux ;
    flux_ph_c_bn       = v.NP(61)	* IkBbNFkB * IKK_flux ;
    flux_ph_c_en       = v.NP(62)	* IkBeNFkB * IKK_flux ;

    % A20 Fluxes (NEW)
    flux_rsu_a20       = v.NP(63) ;
    flux_rsr_a20       = v.NP(64)  * (delayed_nfkbn_a20 ^ v.NP(65)) ;
    flux_rd_a20        = v.NP(67)  * a20t ;
    flux_ps_c_a20      = v.NP(68)  * delayed_a20t ;
    flux_pd_c_a20      = v.NP(70)  * a20  ;

    if( t > v.NP(71))  % disables inducible A20 txn to match exp data
        flux_rsr_a20 = 0;
    end


    % -- IKK Activation Module & TLR4 module
    % LPS binding to the surface
    flux_b_LPS         = v.IP(1)   * LPSo  ;

    % ligand binding
    flux_b_TLR4LPS     = v.IP(5)    * TLR4 * LPS ;
    flux_ub_TLR4LPS    = v.IP(6)    * TLR4LPS ;
    flux_b_TLR4LPSen   = v.IP(7)    * TLR4en * LPSen ;
    flux_ub_TLR4LPSen  = v.IP(8)    * TLR4LPSen ;

    % Activation modified on 0426
    flux_a_MyD88       = v.IP(18)    * (TLR4LPS)^v.IP(19) /((TLR4LPS)^v.IP(19) + (v.IP(20))^v.IP(19))  * MyD88;
    flux_i_MyD88       = v.IP(21)    * MyD88s ;
    flux_a_TRIF        = v.IP(22)    * TLR4LPSen * TRIF;
    flux_i_TRIF        = v.IP(23)   * TRIFs ;

    % shuttling, modified on 0426
    flux_in_LPS        = v.IP(2)   * LPS ;
    flux_out_LPS       = v.IP(3)   * LPSen;
    flux_in_TLR4       = v.IP(12)   * TLR4;
    flux_out_TLR4      = v.IP(13)   * TLR4en;
    flux_in_TLR4LPS    = v.IP(14)   * TLR4LPS;
    flux_out_TLR4LPS   = v.IP(15)   * TLR4LPSen;

    % generation and degradation
    flux_g_TLR4        = v.IP(9);
    flux_d_TLR4        = v.IP(10)   * TLR4;
    flux_d_TLR4en      = v.IP(11)   * TLR4en;
    flux_d_LPSen       = v.IP(4)   * LPSen;
    flux_d_TLR4LPSen   = v.IP(17)   * TLR4LPSen;
    %modified 0515
    flux_d_TLR4LPS     = v.IP(16)   * TLR4LPS;
    %-----------------------IKK module
    % TRAF6
    flux_a_TRAF6_MyD88s= v.IP(24)   * TRAF6 * MyD88s;
    flux_a_TRAF6_TRIFs = v.IP(25)   * TRAF6 * TRIFs;
    flux_i_TRAF6       = v.IP(26)   * TRAF6s;

    % IKKK

    flux_IKKK_on       = v.IP(28)   * IKKK_off;
    flux_IKKK_on_TRAF6s= v.IP(27)   * IKKK_off * TRAF6s;
    flux_IKKK_off      = v.IP(29)   * IKKK;

    % IKK
    flux_IKK_on        = v.IP(31)   * IKK_off;
    flux_IKK_on_IKKK   = v.IP(30)   * IKK_off * IKKK;
    flux_IKK_off       = v.IP(32)   * IKK;
    flux_IKK_off_i     = v.IP(33)   * IKK;
    flux_IKKi          = v.IP(34)   * IKK_i;

    % IRF3  module ------------------------------
    %TBK activation
    flux_a_TBK1_TRIFs  = v.IP(35)   *TBK1*TRIFs;
    flux_i_TBK1        = v.IP(36)   *TBK1s;

    %IRF3 activation
    flux_a_IRF3_TBK1s  = v.IP(37)   *IRF3*TBK1s;
    flux_i_IRF3        = v.IP(38)   *IRF3s;
    flux_i_IRF3n       = v.IP(41)   *IRF3ns;

    %IRF3 shuttling
    flux_out_IRF3      = v.IP(42)   *IRF3n;        
    flux_in_IRF3       = v.IP(43)   *IRF3;
    flux_in_IRF3s      = v.IP(39)   *IRF3s;
    flux_out_IRF3s     = v.IP(40)   *IRF3ns;

    %IRF3 generation and degradation 
    flux_g_IRF3   =  v.IP(44);
    flux_d_IRF3ns =  v.IP(48) * IRF3ns;
    flux_d_IRF3n  =  v.IP(47) * IRF3n;
    flux_d_IRF3   =  v.IP(45) * IRF3;
    flux_d_IRF3s  =  v.IP(46) * IRF3s;

    %TNFR part (NEW)
    flux_pd_m_tnf      = v.IP(53)    * TNF ; % pd_m_tnf 45' half life of exogenous TNF
    flux_syn_tnfrm     = v.IP(54);  % tnfrm synthesis (txn, tsl, localization)
    flux_pd_tnfrm      = v.IP(55)    * tnfrm;% tnfrm --> deg  -- 120' halflife
    flux_a_tnfrm       = v.IP(56)    * tnfrm;% 3tnfrm --> TNFR
    flux_d_TNFR        = v.IP(57)    * TNFR;% TNFR   --> 3tnfrm
    flux_i_TNFR        = v.IP(58)    * TNFR;% TNFR internalization -- 30' halflife
    flux_a_C1_off      = v.IP(59)    * TNFR* TTR;% TNFR + TTR --> C1_off
    flux_d_C1_off      = v.IP(60)    * C1_off;% C1_off --> TNFR + TTR
    flux_a_C1          = v.IP(61)    * C1_off;% C1_off --> C1
    flux_C1_off        = v.IP(62)   * C1;% C1     --> C1_off
    flux_C1_A20        = v.IP(63)   * C1 * a20;% C1     --> C1_off (A20 Mediated)
                                               % flux_d_C1          = v.IP(64)   * C1_off; %ERROR
    flux_d_C1          = v.IP(64)   * C1; %UPDATE % C1     --> TNFR + TTR
    flux_i_C1_off      = v.IP(65)   * C1_off; % C1_off internalization
    flux_i_C1          = v.IP(66)   * C1; % C1 internalization
    flux_a_tnfrmtnf    = v.IP(67)   * tnfrm * TNF; % 3tnfrm + tnf --> TNFRtnf
    flux_a_TNFRtnf     = v.IP(68)   * TNFR * TNF; % TNFR + tnf --> TNFRtnf
    flux_d_TNFRtnf     = v.IP(69)   * TNFRtnf; % TNFRtnf   --> TNFR + tnf
    flux_i_TNFRtnf     = v.IP(70)   * TNFRtnf;  % TNFRtnf internalization -- 30' halflife

    flux_a_C1tnf_off   = v.IP(71)   * TNFRtnf * TTR; % TNFRtnf + TTR --> C1tnf_off
    flux_d_C1tnf_off   = v.IP(72)   * C1tnf_off; % C1tnf_off --> TNFRtnf + TTR
    flux_a_C1tnf       = v.IP(73)   * C1tnf_off; % C1tnf_off --> C1tnf
    flux_C1tnf_off     = v.IP(74)   * C1tnf; % C1tnf     --> C1tnf_off
    flux_C1tnf_A20     = v.IP(75)   * C1tnf * a20;  % C1tnf     --> C1tnf_off (A20 Mediated)
                                                    %flux_d_C1tnf       = v.IP(76)   * C1tnf_off; % ERROR
    flux_d_C1tnf       = v.IP(76)   * C1tnf; % updated  % C1tnf     --> TNFRtnf + TTR
    flux_i_C1tnf       = v.IP(78)   * C1tnf; % C1tnf internalization
    flux_i_C1tnf_off   = v.IP(77)   * C1tnf_off;  % C1tnf_off internalization
    flux_d_tnf_C1_off  = v.IP(79)   * C1tnf_off; % C1tnf_off --> C1_off + tnf
    flux_a_tnf_C1_off  = v.IP(80)   * C1_off * TNF; % C1_off + tnf --> C1tnf_off
    flux_d_tnf_C1      = v.IP(81)   * C1tnf; % C1tnf    --> C1 + tnf
    flux_a_tnf_C1      = v.IP(82)   * C1 * TNF; % C1 + tnf --> C1tnf
    flux_IKKK_on_C1    = v.IP(83)   * IKKK_off * C1; % IKKK_off --> IKKK (C1 mediated)?500
    flux_IKKK_on_C1tnf = v.IP(84)   * IKKK_off * C1tnf;% IKKK_off --> IKKK (C1tnf mediated)

    % CpG-TLR9 (NEW)
    flux_g_tlr9             = v.IP(85);
    flux_d_tlr9             = v.IP(86) * TLR9;
    flux_b_cpg_tlr9         = v.IP(87) * CpG * TLR9;
    flux_d_cpg_tlr9         = v.IP(88) * CpGTLR9;
    flux_a_MyD88_cpgtlr9    = v.IP(89) * CpGTLR9^3/(CpGTLR9^3 + v.IP(90)^3);


    % PIC-TLR3 (NEW)
    flux_g_tlr3             = v.IP(91);
    flux_d_tlr3             = v.IP(92) * TLR3;
    flux_b_pic_tlr3         = v.IP(93) * PIC * TLR3;
    flux_d_pic_tlr3         = v.IP(94) * PICTLR3;
    flux_a_TRIF_pictlr3     = v.IP(95) * TRIF * PICTLR3;



    %% Add Fluxes to component concentrations, then Save

    % IkB Transcrv.IPtion
    delta_IkBat     = delta_IkBat + flux_rsu_a;
    delta_IkBbt     = delta_IkBbt + flux_rsu_b;
    delta_IkBet     = delta_IkBet + flux_rsu_e;
    delta_IkBat     = delta_IkBat + flux_rsr_an;
    delta_IkBbt     = delta_IkBbt + flux_rsr_bn;
    delta_IkBet     = delta_IkBet + flux_rsr_en;
    delta_TNFnas    = delta_TNFnas + flux_rsu_f; %NEW
    delta_TNFnas    = delta_TNFnas + flux_rsr_fn; %NEW
    delta_TNFnas    = delta_TNFnas - flux_rp_f; %NEW
    delta_TNFmRNA   = delta_TNFmRNA + flux_rp_f; %NEW


    % IkB Transcrv.IPt Degradation
    delta_IkBbt     = delta_IkBbt - flux_rd_b;
    delta_IkBet     = delta_IkBet - flux_rd_e;
    delta_IkBat     = delta_IkBat - flux_rd_a;
    delta_TNFmRNA   = delta_TNFmRNA  - flux_rd_f;   %NEW, tnft
                                                    %degradation 

    % IkB Translation
    delta_IkBat     = delta_IkBat - flux_ps_c_a;
    delta_IkBat     = delta_IkBat + flux_ps_c_a;
    delta_IkBa      = delta_IkBa  + flux_ps_c_a;

    delta_IkBbt     = delta_IkBbt - flux_ps_c_b;
    delta_IkBbt     = delta_IkBbt + flux_ps_c_b;
    delta_IkBb      = delta_IkBb  + flux_ps_c_b;

    delta_IkBet     = delta_IkBet - flux_ps_c_e;
    delta_IkBet     = delta_IkBet + flux_ps_c_e;
    delta_IkBe      = delta_IkBe  + flux_ps_c_e;

    delta_TNFpro    = delta_TNFpro + flux_ps_c_f; %NEW 

    % TNF secretion + protein degradation 
    delta_TNFpro     = delta_TNFpro - flux_sc_c_f; %NEW 
    delta_TNF        = delta_TNF    + flux_sc_c_f; %NEW 
    delta_TNFpro     = delta_TNFpro - flux_pd_f; %NEW 

    % IkB:NFkB Shuttling (Free and Bound)
    delta_IkBa      = delta_IkBa  - flux_in_a;
    delta_IkBan     = delta_IkBan + flux_in_a;

    delta_IkBb      = delta_IkBb  - flux_in_b;
    delta_IkBbn     = delta_IkBbn + flux_in_b;

    delta_IkBe      = delta_IkBe  - flux_in_e;
    delta_IkBen     = delta_IkBen + flux_in_e;

    delta_NFkB      = delta_NFkB  - flux_in_n;
    delta_NFkBn     = delta_NFkBn + flux_in_n;

    delta_IkBan     = delta_IkBan - flux_ex_a;
    delta_IkBa      = delta_IkBa  + flux_ex_a;

    delta_IkBbn     = delta_IkBbn - flux_ex_b;
    delta_IkBb      = delta_IkBb  + flux_ex_b;

    delta_IkBen     = delta_IkBen - flux_ex_e;
    delta_IkBe      = delta_IkBe  + flux_ex_e;

    delta_NFkBn     = delta_NFkBn - flux_ex_n;
    delta_NFkB      = delta_NFkB  + flux_ex_n;

    delta_IkBaNFkB  = delta_IkBaNFkB  - flux_in_2an;
    delta_IkBaNFkBn = delta_IkBaNFkBn + flux_in_2an;

    delta_IkBbNFkB  = delta_IkBbNFkB  - flux_in_2bn;
    delta_IkBbNFkBn = delta_IkBbNFkBn + flux_in_2bn;

    delta_IkBeNFkB  = delta_IkBeNFkB  - flux_in_2en;
    delta_IkBeNFkBn = delta_IkBeNFkBn + flux_in_2en;

    delta_IkBaNFkBn = delta_IkBaNFkBn - flux_ex_2an;
    delta_IkBaNFkB  = delta_IkBaNFkB  + flux_ex_2an;

    delta_IkBbNFkBn = delta_IkBbNFkBn - flux_ex_2bn;
    delta_IkBbNFkB  = delta_IkBbNFkB  + flux_ex_2bn;

    delta_IkBeNFkBn = delta_IkBeNFkBn - flux_ex_2en;
    delta_IkBeNFkB  = delta_IkBeNFkB  + flux_ex_2en;

    % IkB:NFkB Association (Cytoplasm and Nucleus)
    delta_IkBa      = delta_IkBa - flux_a_c_an;
    delta_NFkB      = delta_NFkB - flux_a_c_an;
    delta_IkBaNFkB  = delta_IkBaNFkB + flux_a_c_an;

    delta_IkBb      = delta_IkBb - flux_a_c_bn;
    delta_NFkB      = delta_NFkB - flux_a_c_bn;
    delta_IkBbNFkB  = delta_IkBbNFkB + flux_a_c_bn;

    delta_IkBe      = delta_IkBe - flux_a_c_en;
    delta_NFkB      = delta_NFkB - flux_a_c_en;
    delta_IkBeNFkB  = delta_IkBeNFkB + flux_a_c_en;

    delta_IkBan     = delta_IkBan - flux_a_n_an;
    delta_NFkBn     = delta_NFkBn - flux_a_n_an;
    delta_IkBaNFkBn = delta_IkBaNFkBn + flux_a_n_an;

    delta_IkBbn     = delta_IkBbn - flux_a_n_bn;
    delta_NFkBn     = delta_NFkBn - flux_a_n_bn;
    delta_IkBbNFkBn = delta_IkBbNFkBn + flux_a_n_bn;

    delta_IkBen     = delta_IkBen - flux_a_n_en;
    delta_NFkBn     = delta_NFkBn - flux_a_n_en;
    delta_IkBeNFkBn = delta_IkBeNFkBn + flux_a_n_en;

    % IkB:NFkB Dissociation (Cytoplasm and Nucleus)
    delta_IkBaNFkB  = delta_IkBaNFkB - flux_d_c_an;
    delta_IkBa      = delta_IkBa + flux_d_c_an;
    delta_NFkB      = delta_NFkB + flux_d_c_an;

    delta_IkBbNFkB  = delta_IkBbNFkB - flux_d_c_bn;
    delta_IkBb      = delta_IkBb + flux_d_c_bn;
    delta_NFkB      = delta_NFkB + flux_d_c_bn;

    delta_IkBeNFkB  = delta_IkBeNFkB - flux_d_c_en;
    delta_IkBe      = delta_IkBe + flux_d_c_en;
    delta_NFkB      = delta_NFkB + flux_d_c_en;

    delta_IkBaNFkBn = delta_IkBaNFkBn - flux_d_n_an;
    delta_IkBan     = delta_IkBan + flux_d_n_an;
    delta_NFkBn     = delta_NFkBn + flux_d_n_an;

    delta_IkBbNFkBn = delta_IkBbNFkBn - flux_d_n_bn;
    delta_IkBbn     = delta_IkBbn + flux_d_n_bn;
    delta_NFkBn     = delta_NFkBn + flux_d_n_bn;

    delta_IkBeNFkBn = delta_IkBeNFkBn - flux_d_n_en;
    delta_IkBen     = delta_IkBen + flux_d_n_en;
    delta_NFkBn     = delta_NFkBn + flux_d_n_en;

    % Free IkB Degradation (Cytoplasm and Nucleus)
    delta_IkBa      = delta_IkBa  - flux_pd_c_a;
    delta_IkBb      = delta_IkBb  - flux_pd_c_b;
    delta_IkBe      = delta_IkBe  - flux_pd_c_e;
    delta_IkBan     = delta_IkBan - flux_pd_n_a;
    delta_IkBbn     = delta_IkBbn - flux_pd_n_b;
    delta_IkBen     = delta_IkBen - flux_pd_n_e;

    % IkB:NFkB Degradation (Cytoplasm and Nucleus)
    delta_IkBaNFkB  = delta_IkBaNFkB - flux_pd_c_2an;
    delta_NFkB      = delta_NFkB + flux_pd_c_2an;

    delta_IkBbNFkB  = delta_IkBbNFkB - flux_pd_c_2bn;
    delta_NFkB      = delta_NFkB + flux_pd_c_2bn;

    delta_IkBeNFkB  = delta_IkBeNFkB - flux_pd_c_2en;
    delta_NFkB      = delta_NFkB + flux_pd_c_2en;

    delta_IkBaNFkBn = delta_IkBaNFkBn - flux_pd_n_2an;
    delta_NFkBn     = delta_NFkBn + flux_pd_n_2an;

    delta_IkBbNFkBn = delta_IkBbNFkBn - flux_pd_n_2bn;
    delta_NFkBn     = delta_NFkBn + flux_pd_n_2bn;

    delta_IkBeNFkBn = delta_IkBeNFkBn - flux_pd_n_2en;
    delta_NFkBn     = delta_NFkBn + flux_pd_n_2en;

    % IKK Mediated IkB Degradation
    delta_IkBa      = delta_IkBa - flux_ph_c_a;
    delta_IkBb      = delta_IkBb - flux_ph_c_b;
    delta_IkBe      = delta_IkBe - flux_ph_c_e;

    % IKK Mediated IkB:NFkB Degradation
    delta_IkBaNFkB  = delta_IkBaNFkB - flux_ph_c_an;
    delta_NFkB      = delta_NFkB + flux_ph_c_an;

    delta_IkBbNFkB  = delta_IkBbNFkB - flux_ph_c_bn;
    delta_NFkB      = delta_NFkB + flux_ph_c_bn;

    delta_IkBeNFkB  = delta_IkBeNFkB - flux_ph_c_en;
    delta_NFkB      = delta_NFkB + flux_ph_c_en;


    % Upstream Pathway  [LPS] --> IKK activation

    % Ligand binding
    delta_LPSo      = delta_LPSo - flux_b_LPS;
    delta_LPS       = delta_LPS  + flux_b_LPS;

    delta_LPS       = delta_LPS       -   flux_b_TLR4LPS;
    delta_TLR4       = delta_TLR4       -   flux_b_TLR4LPS;
    delta_TLR4LPS   = delta_TLR4LPS   +   flux_b_TLR4LPS;
    delta_LPS       = delta_LPS       +   flux_ub_TLR4LPS;
    delta_TLR4       = delta_TLR4     +   flux_ub_TLR4LPS;
    delta_TLR4LPS   = delta_TLR4LPS   -   flux_ub_TLR4LPS;

    delta_LPSen     = delta_LPSen     -   flux_b_TLR4LPSen;
    delta_TLR4en     = delta_TLR4en     -   flux_b_TLR4LPSen;
    delta_TLR4LPSen = delta_TLR4LPSen +   flux_b_TLR4LPSen;
    delta_LPSen     = delta_LPSen     +   flux_ub_TLR4LPSen;
    delta_TLR4en     = delta_TLR4en   +   flux_ub_TLR4LPSen;
    delta_TLR4LPSen = delta_TLR4LPSen -   flux_ub_TLR4LPSen;

    % Activation
    delta_MyD88     = delta_MyD88     -   flux_a_MyD88;
    delta_MyD88s    = delta_MyD88s    +   flux_a_MyD88;
    delta_MyD88     = delta_MyD88     +   flux_i_MyD88;
    delta_MyD88s    = delta_MyD88s    -   flux_i_MyD88;

    delta_TRIF      = delta_TRIF      -   flux_a_TRIF;
    delta_TRIFs     = delta_TRIFs     +   flux_a_TRIF;
    delta_TRIF      = delta_TRIF      +   flux_i_TRIF;
    delta_TRIFs     = delta_TRIFs     -   flux_i_TRIF;

    % shuttling
    delta_LPS       = delta_LPS       -   flux_in_LPS;
    delta_LPSen     = delta_LPSen     +   flux_in_LPS;
    delta_LPSo       = delta_LPSo     +   flux_out_LPS;
    delta_LPSen     = delta_LPSen     -   flux_out_LPS;

    delta_TLR4      = delta_TLR4      -   flux_in_TLR4;
    delta_TLR4en    = delta_TLR4en    +   flux_in_TLR4;
    delta_TLR4      = delta_TLR4      +   flux_out_TLR4;
    delta_TLR4en    = delta_TLR4en    -   flux_out_TLR4;

    delta_TLR4LPS      = delta_TLR4LPS      -   flux_in_TLR4LPS;
    delta_TLR4LPSen    = delta_TLR4LPSen    +   flux_in_TLR4LPS;
    delta_TLR4LPS      = delta_TLR4LPS      +   flux_out_TLR4LPS;
    delta_TLR4LPSen     = delta_TLR4LPSen   -   flux_out_TLR4LPS;


    % generation and degradation
    delta_TLR4      = delta_TLR4 + flux_g_TLR4;
    delta_TLR4      = delta_TLR4 - flux_d_TLR4;
    delta_TLR4en     = delta_TLR4en - flux_d_TLR4en;
    delta_LPSen      = delta_LPSen - flux_d_LPSen;
    delta_TLR4LPSen      = delta_TLR4LPSen - flux_d_TLR4LPSen;
    %-- modified 0515
    delta_TLR4LPS      = delta_TLR4LPS - flux_d_TLR4LPS;


    % TRAF6
    delta_TRAF6       = delta_TRAF6  - flux_a_TRAF6_MyD88s;
    delta_TRAF6s      = delta_TRAF6s + flux_a_TRAF6_MyD88s;
    delta_TRAF6       = delta_TRAF6  - flux_a_TRAF6_TRIFs;
    delta_TRAF6s      = delta_TRAF6s + flux_a_TRAF6_TRIFs;
    delta_TRAF6       = delta_TRAF6  + flux_i_TRAF6;
    delta_TRAF6s      = delta_TRAF6s - flux_i_TRAF6;

    % IKKK Regulation
    delta_IKKK      = delta_IKKK     + flux_IKKK_on;
    delta_IKKK_off  = delta_IKKK_off - flux_IKKK_on;

    delta_IKKK      = delta_IKKK     + flux_IKKK_on_TRAF6s;
    delta_IKKK_off  = delta_IKKK_off - flux_IKKK_on_TRAF6s;

    delta_IKKK_off  = delta_IKKK_off + flux_IKKK_off;
    delta_IKKK      = delta_IKKK     - flux_IKKK_off;

    % IKK Regulation
    delta_IKK       = delta_IKK     + flux_IKK_on;
    delta_IKK_off   = delta_IKK_off - flux_IKK_on;

    delta_IKK       = delta_IKK     + flux_IKK_on_IKKK;
    delta_IKK_off   = delta_IKK_off - flux_IKK_on_IKKK;

    delta_IKK       = delta_IKK     - flux_IKK_off;
    delta_IKK_off   = delta_IKK_off + flux_IKK_off;

    delta_IKK       = delta_IKK     - flux_IKK_off_i;
    delta_IKK_i     = delta_IKK_i   + flux_IKK_off_i;

    delta_IKK_off   = delta_IKK_off + flux_IKKi;
    delta_IKK_i     = delta_IKK_i   - flux_IKKi;

    % IRF3 module
    % activation & inactivation
    delta_TBK1      = delta_TBK1 + flux_i_TBK1;
    delta_TBK1s     = delta_TBK1s - flux_i_TBK1;
    delta_TBK1      = delta_TBK1 - flux_a_TBK1_TRIFs;
    delta_TBK1s     = delta_TBK1s + flux_a_TBK1_TRIFs;
    delta_IRF3      = delta_IRF3 + flux_i_IRF3;
    delta_IRF3s     = delta_IRF3s - flux_i_IRF3;
    delta_IRF3      = delta_IRF3 - flux_a_IRF3_TBK1s;
    delta_IRF3s     = delta_IRF3s + flux_a_IRF3_TBK1s;
    delta_IRF3n      = delta_IRF3n + flux_i_IRF3n;
    delta_IRF3ns     = delta_IRF3ns - flux_i_IRF3n;

    % shuttling.
    delta_IRF3n      = delta_IRF3n + flux_in_IRF3;
    delta_IRF3       = delta_IRF3  - flux_in_IRF3;
    delta_IRF3n      = delta_IRF3n - flux_out_IRF3;
    delta_IRF3       = delta_IRF3  + flux_out_IRF3;
    delta_IRF3ns      = delta_IRF3ns + flux_in_IRF3s;
    delta_IRF3s       = delta_IRF3s  - flux_in_IRF3s;
    delta_IRF3ns      = delta_IRF3ns - flux_out_IRF3s;
    delta_IRF3s       = delta_IRF3s  + flux_out_IRF3s;


    % generation and degradation
    delta_IRF3        = delta_IRF3   + flux_g_IRF3;
    delta_IRF3ns      = delta_IRF3ns - flux_d_IRF3ns;
    delta_IRF3n       = delta_IRF3n  - flux_d_IRF3n;
    delta_IRF3s       = delta_IRF3s  - flux_d_IRF3s;
    delta_IRF3        = delta_IRF3   - flux_d_IRF3;


    %% TNR module (NEW) + A20 (NEW)
    % A20 Expression
    % Transcrv.IPtion
    delta_a20t      = delta_a20t + flux_rsu_a20;
    delta_a20t      = delta_a20t + flux_rsr_a20;
    delta_a20t      = delta_a20t - flux_rd_a20;

    % Translation
    delta_a20       = delta_a20  + flux_ps_c_a20;

    % Degradation
    delta_a20       = delta_a20 - flux_pd_c_a20;

    % Upstream Pathway  [TNF] --> IKK activation

    % TNF-independent Activation

    % tnfrm trimerization 3tnfrm<--> TNFR
    delta_tnfrm     = delta_tnfrm- 3*flux_a_tnfrm;
    delta_TNFR      = delta_TNFR +   flux_a_tnfrm;

    delta_tnfrm     = delta_tnfrm+ 3*flux_d_TNFR;
    delta_TNFR      = delta_TNFR -   flux_d_TNFR;

    % TNF Receptor Metabolism
    delta_TNFR     = delta_TNFR     - flux_i_TNFR;
    delta_TNFRtnf  = delta_TNFRtnf  - flux_i_TNFRtnf;

    delta_tnfrm    = delta_tnfrm    - flux_pd_tnfrm;
    delta_tnfrm    = delta_tnfrm    + flux_syn_tnfrm;

    % Complex I Generation  TTR + TNFR <--> C1_off, C1
    delta_C1_off    = delta_C1_off  + flux_a_C1_off;
    delta_TTR       = delta_TTR     - flux_a_C1_off;
    delta_TNFR      = delta_TNFR    - flux_a_C1_off;

    delta_C1_off    = delta_C1_off  - flux_d_C1_off;
    delta_TTR       = delta_TTR     + flux_d_C1_off;
    delta_TNFR      = delta_TNFR    + flux_d_C1_off;

    delta_C1        = delta_C1      - flux_d_C1;
    delta_TTR       = delta_TTR     + flux_d_C1;
    delta_TNFR      = delta_TNFR    + flux_d_C1;

    % C1_off  <--> C1
    delta_C1_off    = delta_C1_off  - flux_a_C1;
    delta_C1        = delta_C1      + flux_a_C1;

    delta_C1_off    = delta_C1_off  + flux_C1_A20;
    delta_C1        = delta_C1      - flux_C1_A20;

    delta_C1_off    = delta_C1_off  + flux_C1_off;
    delta_C1        = delta_C1      - flux_C1_off;

    % C1 and C1_off internalization (treated as deg)
    delta_C1        = delta_C1      - flux_i_C1;
    delta_C1_off    = delta_C1_off  - flux_i_C1_off;

    % TNF-Dependent Activation

    % TNF degradation
    delta_TNF       = delta_TNF     - flux_pd_m_tnf;

    % TNF:TNFR binding
    delta_TNFRtnf      = delta_TNFRtnf    +   flux_a_tnfrmtnf;%UPDATE
                                                              %         delta_TNFR      = delta_TNFR    +   flux_a_tnfrmtnf; %ERROR
    delta_tnfrm     = delta_tnfrm   - 3*flux_a_tnfrmtnf;
    delta_TNF       = delta_TNF     -   flux_a_tnfrmtnf;

    delta_TNFR      = delta_TNFR    - flux_a_TNFRtnf;
    delta_TNF       = delta_TNF     - flux_a_TNFRtnf;
    delta_TNFRtnf   = delta_TNFRtnf + flux_a_TNFRtnf;

    delta_TNFR      = delta_TNFR    + flux_d_TNFRtnf;
    delta_TNF       = delta_TNF     + flux_d_TNFRtnf;
    delta_TNFRtnf   = delta_TNFRtnf - flux_d_TNFRtnf;

    % Complex I Generation  TTR + TNFRtnf <--> C1tnf_off, C1tnf
    delta_C1tnf_off = delta_C1tnf_off   + flux_a_C1tnf_off;
    delta_TTR       = delta_TTR         - flux_a_C1tnf_off;
    delta_TNFRtnf   = delta_TNFRtnf     - flux_a_C1tnf_off;

    delta_C1tnf_off = delta_C1tnf_off   - flux_d_C1tnf_off;
    delta_TTR       = delta_TTR         + flux_d_C1tnf_off;
    delta_TNFRtnf   = delta_TNFRtnf     + flux_d_C1tnf_off;

    delta_C1tnf     = delta_C1tnf       - flux_d_C1tnf;
    delta_TTR       = delta_TTR         + flux_d_C1tnf;
    delta_TNFRtnf   = delta_TNFRtnf     + flux_d_C1tnf;

    % C1tnf_off <--> C1tnf
    delta_C1tnf_off = delta_C1tnf_off   - flux_a_C1tnf;
    delta_C1tnf     = delta_C1tnf       + flux_a_C1tnf;

    delta_C1tnf_off = delta_C1tnf_off   + flux_C1tnf_off;
    delta_C1tnf     = delta_C1tnf       - flux_C1tnf_off;

    delta_C1tnf_off = delta_C1tnf_off   + flux_C1tnf_A20;
    delta_C1tnf     = delta_C1tnf       - flux_C1tnf_A20;

    % C1tnf and C1tnf_off internalization (treated as deg)
    delta_C1tnf_off = delta_C1tnf_off   - flux_i_C1tnf_off;
    delta_C1tnf     = delta_C1tnf       - flux_i_C1tnf;

    % C1tnf_off,C1tnf <--> C1_off,C1 + TNF
    delta_C1tnf_off = delta_C1tnf_off   - flux_d_tnf_C1_off;
    delta_C1_off    = delta_C1_off      + flux_d_tnf_C1_off;
    delta_TNF       = delta_TNF         + flux_d_tnf_C1_off;

    delta_C1tnf_off = delta_C1tnf_off   + flux_a_tnf_C1_off;
    delta_C1_off    = delta_C1_off      - flux_a_tnf_C1_off;
    delta_TNF       = delta_TNF         - flux_a_tnf_C1_off;

    delta_C1tnf     = delta_C1tnf       - flux_d_tnf_C1;
    delta_C1        = delta_C1          + flux_d_tnf_C1;
    delta_TNF       = delta_TNF         + flux_d_tnf_C1;

    delta_C1tnf     = delta_C1tnf       + flux_a_tnf_C1;
    delta_C1        = delta_C1          - flux_a_tnf_C1;
    delta_TNF       = delta_TNF         - flux_a_tnf_C1;

    % IKKK Regulation

    delta_IKKK      = delta_IKKK     + flux_IKKK_on_C1;
    delta_IKKK_off  = delta_IKKK_off - flux_IKKK_on_C1;

    delta_IKKK      = delta_IKKK     + flux_IKKK_on_C1tnf;
    delta_IKKK_off  = delta_IKKK_off - flux_IKKK_on_C1tnf;

    %% TLR3 (NEW) + TLR9 (NEW)
    % CpG-TLR9 (NEW)
    delta_TLR9      = delta_TLR9 + flux_g_tlr9;
    delta_TLR9      = delta_TLR9 - flux_d_tlr9;

    delta_TLR9      = delta_TLR9    - flux_b_cpg_tlr9;
    delta_CpG       = delta_CpG     - flux_b_cpg_tlr9;
    delta_CpGTLR9   = delta_CpGTLR9 + flux_b_cpg_tlr9;

    delta_TLR9      = delta_TLR9    + flux_d_cpg_tlr9;
    delta_CpG       = delta_CpG     + flux_d_cpg_tlr9;
    delta_CpGTLR9   = delta_CpGTLR9 - flux_d_cpg_tlr9;

    delta_MyD88s    = delta_MyD88s  + flux_a_MyD88_cpgtlr9;
    delta_MyD88     = delta_MyD88   - flux_a_MyD88_cpgtlr9;

    % PIC-TLR3 (NEW)
    delta_TLR3      = delta_TLR3 + flux_g_tlr3;
    delta_TLR3      = delta_TLR3 - flux_d_tlr3;

    delta_TLR3      = delta_TLR3    - flux_b_pic_tlr3;
    delta_PIC       = delta_PIC     - flux_b_pic_tlr3;
    delta_PICTLR3   = delta_PICTLR3 + flux_b_pic_tlr3;

    delta_TLR3      = delta_TLR3    + flux_d_pic_tlr3;
    delta_PIC       = delta_PIC     + flux_d_pic_tlr3;
    delta_PICTLR3   = delta_PICTLR3 - flux_d_pic_tlr3;

    delta_TRIFs    = delta_TRIFs  + flux_a_TRIF_pictlr3;
    delta_TRIF     = delta_TRIF   - flux_a_TRIF_pictlr3;

    %==========================================================
    %% Save concentrations for next time step
    delta(1,1)      = delta_IkBa;
    delta(2,1)      = delta_IkBan;
    delta(3,1)      = delta_IkBaNFkB;
    delta(4,1)      = delta_IkBaNFkBn;
    delta(5,1)      = delta_IkBat;
    delta(6,1)      = delta_IkBb;
    delta(7,1)      = delta_IkBbn;
    delta(8,1)      = delta_IkBbNFkB;
    delta(9,1)      = delta_IkBbNFkBn;
    delta(10,1)     = delta_IkBbt;
    delta(11,1)     = delta_IkBe;
    delta(12,1)     = delta_IkBen;
    delta(13,1)     = delta_IkBeNFkB;
    delta(14,1)     = delta_IkBeNFkBn;
    delta(15,1)     = delta_IkBet;
    delta(16,1)     = delta_NFkB;
    delta(17,1)     = delta_NFkBn;
    %
    delta(18,1)     = delta_LPSo;        
    delta(19,1)     = delta_LPS;
    delta(20,1)     = delta_LPSen;        
    delta(21,1)     = delta_TLR4;
    delta(22,1)     = delta_TLR4en;        
    delta(23,1)     = delta_TLR4LPS;
    delta(24,1)     = delta_TLR4LPSen;        
    delta(25,1)     = delta_MyD88;
    delta(26,1)     = delta_MyD88s;
    delta(27,1)     = delta_TRIF;
    delta(28,1)     = delta_TRIFs;
    delta(29,1)     = delta_TRAF6;
    delta(30,1)     = delta_TRAF6s;
    delta(31,1)     = delta_IKKK_off;        
    delta(32,1)     = delta_IKKK;
    delta(33,1)     = delta_IKK_off;
    delta(34,1)     = delta_IKK;
    delta(35,1)     = delta_IKK_i;

    %IRF module
    delta(36,1)     = delta_TBK1;
    delta(37,1)     = delta_TBK1s;
    delta(38,1)     = delta_IRF3;
    delta(39,1)     = delta_IRF3s;        
    delta(40,1)     = delta_IRF3n;
    delta(41,1)     = delta_IRF3ns;

    % TNF component (NEW)
    delta(42,1)     = delta_TNFnas  ;
    delta(43,1)     = delta_TNFmRNA ;
    delta(44,1)     = delta_TNFpro  ;
    delta(45,1)     = delta_TNF  ;        

    % TNFR part (NEW, total 8)
    delta(46,1)     = delta_tnfrm;  
    delta(47,1)     = delta_TNFR;         
    delta(48,1)     = delta_TNFRtnf;         
    delta(49,1)     = delta_C1;         
    delta(50,1)     = delta_C1_off  ;         
    delta(51,1)     = delta_C1tnf  ;  
    delta(52,1)     = delta_C1tnf_off   ;  
    delta(53,1)     = delta_TTR ;  

    % A20 (NEW) 
    delta(54,1)     = delta_a20 ;  
    delta(55,1)     = delta_a20t ;  

    % CpG-TLR9 (NEW)
    delta(56,1)     = delta_CpG  ;       
    delta(57,1)     = delta_TLR9 ;
    delta(58,1)     = delta_CpGTLR9;

    % PIC-TLR3 (NEW)
    delta(59,1)     = delta_PIC;
    delta(60,1)     = delta_TLR3;
    delta(61,1)     = delta_PICTLR3;
end