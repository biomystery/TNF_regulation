function v=nfkbRunPar(v)
%% Intro
% Input: 1. v.LPS_DOSE
%        2. v.IP (wt, MyD88ko, TrifKo) 
% Require: nfkbBasal.m,nfkbSim.m
% Output: v.OUTPUT

%% init
v.START_TIME   = -4000; % min of equilibration phase (phase 1)

%   --- you don't need to change START_TIME
% v.SIM_TIME     = 120;   % min of stimulation phase (phase 2)

% LPS stimulation variables
%    LPS stimulation is accomplished by adding a concentration of LPS
%    ligand at the start of phase 2.
v.LPS_TIME = v.SIM_TIME+1; % min of LPS stimulation. Chronic is SIM_TIME+1
%v.LPS_DOSE = 1e-5;      % missing explanation, suppose 0.05
% uM=100ng/ml, ie. lps molecular weight =
% 2 KDa

% Set Initial Concentrations for Equilibrium Phase
%         v.INITVALUES = getInit();
v.IKK_TOTAL  = sum(v.INITVALUES{1}([33:35]));

% Set Rate Constants (model parameters)
%[v.NP,v.IP] = getRateParams();


%% Call the SIM Function
v = nfkbBasal(v);
v = nfkbSim(v);
end