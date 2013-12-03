function simData=getSimData(id)
%% construct and passing parameters
v=struct;
[v.NP,v.IP] = getRateParams();%wt
v.D_FLAG=1; % display flag
v.P_FLAG=0; % plot flag
v.L_FLAG=0; % Legend turn on or off

% otal time
v.STIMULI  = id.stimuli; % 'TNF' or 'LPS'

v.SIM_TIME     = 360; % min of stimulation phase (phase 2)v.GENOTYPE = id.genotype;

% Stimuli
if strcmp(id.stimuli,'LPS')
    v.DOSE = id.dose*2;
elseif strcmp(id.stimuli,'TNF')
    v.DOSE = id.dose*(1.96e-4);      % 1.96e-4uM= 1ng/mL TNF;
elseif strcmp(id.stimuli,'CpG')
    v.DOSE = id.dose; % uM. start at 0.1 uM
elseif strcmp(id.stimuli,'PIC')
    v.DOSE = id.dose; % uM. start at 0.1 uM
end

v.INITVALUES = getInit();

% mutant
switch id.genotype
    case 'mko'
        v.INITVALUES{1}(21)=0; %myd88
    case 'tko'
        v.INITVALUES{1}(23)=0;
    case 'wt'
        v.IP(24:25) = v.IP(52:53);
end

% Combine parameters for LPS, TLR4 binding.
v.IP(9) = v.IP(7);
v.IP(10)= v.IP(8);

% run
v  = nfkbRunPar(v) ;

%% return
if numel(id.output)==1
    index =find(strcmp(v.INITVALUES{2}(:),id.output));
    simData =v.OUTPUT(:,index)';
elseif numel(id.output)>1
    for i=1:numel(id.output)
        index(i) =find(strcmp(v.INITVALUES{2}(:),id.output{i}));
    end
    simData =v.OUTPUT(:,index)';
else
    error('func:getSim','Wrong output id nubmer');
end
end