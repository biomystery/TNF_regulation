%% basic setting for fig5c
id.genotype = 'wt';
id.flag_noTnfFeedback = true;%true or false
id.output ={'MyD88s','TRIFs','IKKK','IKK','NFkBn','TNFt','TNF'};

%% run the 3 stimulis one by one 
stimulis.name= {'LPS','CpG','PIC'}
stimulis.dose= [10,100,50] % 10ng/ml, 100nM, 50 ug/ml 
colors = {'k'};

for j = 1:3
    id.stimuli = stimuli.name{j}; %LPS, CpG, PIC,TNF. 
    id.dose = stimuli.dose(j); %'1','100' 
    
    sim{j} = getSimData(id);    
    % show the result 
    figure('position',[680   419   217*4   559/3*4])
    for i=1:size(id.output,2)
        subplot(4,3,i);
        plot(sim{1}(i,:),colors{mod(i,1)+1},'linewidth',2);
        xlim([0 360]);set(gca,'xtick',0:60:360);
        title(id.output{i})
    end
end
%% export the data 
nfkbn     = [sim{1}(1,:);sim{2}(1,:);sim{3}(1,:)];
tnfNas    = [sim{1}(2,:);sim{2}(2,:);sim{3}(2,:)];
tnfmRNA   = [sim{1}(3,:);sim{2}(3,:);sim{3}(3,:)];
tnfPro    = [sim{1}(4,:);sim{2}(4,:);sim{3}(4,:)];
tnfSec    = cumsum([sim{1}(4,:)*.18;sim{2}(4,:)*.18;sim{3}(4,:)*.18/1.5]);
IKK    = [sim{1}(5,:);sim{2}(5,:);sim{3}(5,:)];

t = 0:240;
csvwrite('./simData/nfkb_sim.csv',[t;nfkbn]')
csvwrite('./simData/nascent_sim.csv',[t;tnfNas]')
csvwrite('./simData/mRNA_sim.csv',[t;tnfmRNA]')
csvwrite('./simData/prot_sim.csv',[t;tnfPro]')
csvwrite('./simData/sec_sim.csv',[t;tnfSec]')


%% R code to draw the result. 
!R CMD BATCH Fig5c_LPS_CpG_pIC.R
