close all
%% basic setting for fig5c
id.genotype = 'wt';
id.flag_noTnfFeedback = true;%true or false
                             %id.output ={'PIC','CpG','PICTLR3','CpGTLR9','MyD88s',...
                             %    'TRIFs','IKKK','IKK','NFkBn','TNFmRNA','TNF'};
id.output ={'TNFmRNA','TNFpro','NFkBn','IKK','MyD88s','TRIFs','CpG','CpGTLR9'};
id.DT = 1;
id.sim_time = 240;
plot_flag = 1 ;
%% run the 3 stimulis one by one 
stimuli.name= {'LPS','CpG','PIC','TNF'};
stimuli.dose= [100,100,50,1]; % 10ng/ml, 100nM, 50 ug/ml 
colors = {'k'};


for j = 1:3
    id.stimuli = stimuli.name{j}; %LPS, CpG, PIC,TNF. 
    id.dose = stimuli.dose(j); %'1','100' 
    
    sim{j} = getSimData(id);    
    % show the result 
    if plot_flag 
        figure('position',[680   419   217*4   559/3*4])
        for i=1:size(id.output,2)
            subplot(4,3,i);
            plot(sim{j}(i,:),colors{mod(i,1)+1},'linewidth',2);
            xlim([0 id.sim_time]);set(gca,'xtick',0:60:360);
            title(id.output{i})
        end
    end
end

%%
[~,~,tp] = getRateParams(); % LPS , wt 
k_sec = tp(9); %pic condtion, sec rate 
tnfmRNA   = [sim{1}(1,:);sim{2}(1,:);sim{3}(1,:)];
tnfPro    = [sim{1}(2,:);sim{2}(2,:);sim{3}(2,:)];
tnfSec    = cumsum([sim{1}(2,:)*tp(9);sim{2}(2,:)*tp(9)/2.5;sim{3}(2,:)*tp(9)]');


t = 0:id.sim_time;
csvwrite('./simData/mRNA_sim.csv',[t;tnfmRNA]')
csvwrite('./simData/sec_sim.csv',[t;tnfSec']')

% end
!R CMD BATCH Fig5c_LPS_CpG_pIC.R
!rm *.Rout
%!rm *.RD* 
%!rm *.Rhi*