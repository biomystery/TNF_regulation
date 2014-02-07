function simData=getSimData(id)
%% construct and passing parameters
v=struct;
[v.NP,v.IP,v.TP] = getRateParams();%wt
v.D_FLAG=0; % display flag
v.P_FLAG=0; % plot flag
v.L_FLAG=0; % Legend turn on or off
v.DT = id.DT;

% total time
v.SIM_TIME     = id.sim_time; % min of stimulation phase (phase
                              % 2)v.GENOTYPE = id.genotype;
v.STIMULI  = id.stimuli; 

% Stimuli
if strcmp(id.stimuli,'LPS')
    v.DOSE = id.dose*v.IP(52);
    %    v.TP(4) = 0.05;% EC50 nfkb activate tnf
    %    v.TP(2) = 0.05 ;% EC50 nfkb activate tnf            
elseif strcmp(id.stimuli,'TNF')
    v.DOSE = id.dose*(1.96e-4);      % 1.96e-4uM= 1ng/mL TNF;

elseif strcmp(id.stimuli,'CpG') %tko
    v.DOSE = id.dose*(1.96e-4)*100; % no transmition yet
    v.TP(2) = 0.03 ;% EC50 nfkb activate tnf        
    v.TP(4) = 0.02;% EC50 nfkb activate tnf
    v.TP(5) = .07;         % NEW mRNA stability
    v.TP(6) = v.TP(6)/1.5; % NEW process rate 
    v.TP(9) = v.TP(9)/5; % NEW sec rate 
    v.TP(7) = v.TP(7)/1.5; % NEW tl rate     
    v.IP(90) = v.IP(90); % EC50 for MyD88 activation     
    v.IP(87)   = .6;     % CpG +  TLR9 -> CpGTLR9, affinity : 186 nM ± 35 nM
    v.IP(88)   = 2.7; % CpGTLR9 -> CpG +  TLR9 1 × 104 M-1 s-1, 6 mins half-life
    
    
elseif strcmp(id.stimuli,'PIC') % mko 
    v.DOSE = id.dose*(1.96e-4); % no transmistion yet
    v.TP(5) = .07;         % No stabilization for PIC 
    v.TP(6) = v.TP(6)/4.2; % NEW process rate    
    v.TP(4) = v.TP(4) ;% 10 fold higher threshold. 
end


v.INITVALUES = getInit();

switch id.genotype
  case 'mko' %PIC

    v.INITVALUES{1}(25)=0; %myd88
    v.TP(6) = v.TP(6)/4.2; % NEW process rate 4.2
    v.TP(4) = v.TP(4) * 1.5;% 10 fold higher threshold.     
  case 'tko' % CpG

    v.INITVALUES{1}(27)=0;
    v.TP(5) = .07;         % NEW mRNA stability
    v.TP(6) = v.TP(6)/1.5; % NEW process rate 
    v.TP(9) = v.TP(9)/5; % NEW sec rate 
    v.TP(7) = v.TP(7)/1.5; % NEW tl rate     
  case 'wt'
    v.IP(24:25) = v.IP(49:50);
end
% TNF feedback: controlled by TNF receptor  
v.flag_noTnfFeedback = id.flag_noTnfFeedback;
if id.flag_noTnfFeedback
    v.IP(54) = 0; 
end


% Run
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