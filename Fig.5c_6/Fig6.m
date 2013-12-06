%% basic setting for fig5c
id.genotype = 'wt';
%id.output ={'PIC','CpG','PICTLR3','CpGTLR9','MyD88s',...
%    'TRIFs','IKKK','IKK','NFkBn','TNFmRNA','TNF'};
id.output ={'NFkBn'};
id.DT = 1;
id.sim_time = 480;


stimuli.name= {'LPS','CpG','PIC','TNF'};
stimuli.dose= [100,100,50,1]; % 10ng/ml, 100nM, 50 ug/ml 
%% run the 3 stimulis one by one without feedback
id.flag_noTnfFeedback = true;%true or false

for j = 1:3
    id.stimuli = stimuli.name{j}; %LPS, CpG, PIC,TNF. 
    id.dose = stimuli.dose(j); %'1','100' 

    sim{j} = getSimData(id);    

end


%% run the 3 stimulis one by one with feedback
id.flag_noTnfFeedback = false;%true or false

for j = 1:3
    id.stimuli = stimuli.name{j}; %LPS, CpG, PIC,TNF. 
    id.dose = stimuli.dose(j); %'1','100' 

    sim{j+3} = getSimData(id);    

end


%%

nfkbn_nofeedback  = [sim{1}(1,:)/sim{1}(1,1);sim{2}(1,:)/sim{2}(1,1);sim{3}(1,:)/sim{3}(1,1)];
nfkbn_feedback    = [sim{4}(1,:)/sim{4}(1,1);sim{5}(1,:)/sim{5}(1,1);sim{6}(1,:)/sim{6}(1,1)];

t = 0:id.sim_time;
csvwrite('./simData/nfkb_sim_nofeedback.csv',[t;nfkbn_nofeedback]')
csvwrite('./simData/nfkb_sim_feedback.csv',[t;nfkbn_feedback]')

%% end
!R CMD BATCH Fig6.R