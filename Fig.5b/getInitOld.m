function values = getInit()
values     = cell(2,1); % initial concentrations and names

values{1}  = zeros(45,1); % all species start at zero by default
values{2}  = {...
    'IkBa'     ,'IkBan' ,'IkBaNFkB','IkBaNFkBn','IkBat',  ... % 1-5
    'IkBb'     ,'IkBbn' ,'IkBbNFkB','IkBbNFkBn','IkBbt',  ... % 6-10
    'IkBe'     ,'IkBen' ,'IkBeNFkB','IkBeNFkBn','IkBet',  ... % 11-15
    'NFkB'     ,'NFkBn' , 'LPS'    ,'TLR4'     ,'TLR4LPS',... % 16-20
    'MyD88'    ,'MyD88s','TRIF'    ,'TRIFs'    ,'LPSen',  ... % 21-25
    'TLR4LPSen', 'IKKK' , 'IKKK_off','IKK'     ,'IKK_off',... % 26-30,MyD88s means M6
    'IKK_i'    ,'TRAF6' , 'TRAF6s' ,'TLR4en'   ,'TBK1',   ... % 31-35
    'TBK1s'    ,'IRF3'  ,'IRF3n'    ,'IRF3s'   ,'IRF3ns', ... % 36-40
    'LPSo'     ,'TNFnas', 'TNFmRNA',  'TNFpro' ,'TNFsec'  ... % 41-45 NEW
    };

% Set initial values for non-zero species
%   These cycle between different forms (bound/unbound), but
%     the sum always remains at these constant levels
values{1}(17)   = 0.125;  % nfkbn
values{1}(21)   = .1 *100 ;    % MyD88_off
values{1}(23)   = 0.1*100;    % TRIF_off
values{1}(28)   = 0.1*100;    % IKKK_off
values{1}(30)   = 0.1*100;    % IKK_off
values{1}(32)   = 0.1*100;   % TRAF_off
values{1}(35)   = 0.1*100;    % TBK1
values{1}(37)   = 0.1;    % IRF3

% Notes:ss
% - [TLR4] is not held constant, but should
%     equilibrate to  uM (measured?????)
% - Generally, you should start with all IKK in the IKK_off form
end
