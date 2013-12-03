function values = getInit()
values     = cell(2,1); % initial concentrations and names

values{1}  = zeros(44,1); % all species start at zero by default
values{2}  = {...
    'IkBa','IkBan','IkBaNFkB','IkBaNFkBn','IkBat', ... % 1-5
    'IkBb','IkBbn','IkBbNFkB','IkBbNFkBn','IkBbt', ... % 6-10
    'IkBe','IkBen','IkBeNFkB','IkBeNFkBn','IkBet', ... % 11-15
    'NFkB','NFkBn',                                ... % 16-17
    'LPSo','LPS','LPSen','TLR4','TLR4en',          ... % 18-22
    'TLR4LPS','TLR4LPSen','MyD88','MyD88s','TRIF', ... % 23-27,MyD88s means M6
    'TRIFs','TRAF6','TRAF6s','IKKK_off','IKKK',    ... % 28-32
    'IKK_off','IKK','IKK_i','TBK1','TBK1s'         ... % 33-37
    'IRF3','IRF3s','IRF3n','IRF3ns',              ... % 38-41
    'TNFnas', 'TNFmRNA',  'TNFpro'        ...% 42-45
             };

% Set initial values for non-zero species
%   These cycle between different forms (bound/unbound), but
%     the sum always remains at these constant levels
values{1}(17)   = 0.125;  % nfkbn
values{1}(25)   = .1 ;    % MyD88_off
values{1}(27)   = 0.1;    % TRIF_off
values{1}(29)   = 0.1;   % TRAF_off
values{1}(31)   = 0.1;    % IKKK_off
values{1}(33)   = 0.1;    % IKK_off
values{1}(36)   = 0.1;    % TBK1

% Notes:ss
% - [TLR4] is not held constant, but should
%     equilibrate to  uM (measured?????)
% - Generally, you should start with all IKK in the IKK_off form
end
