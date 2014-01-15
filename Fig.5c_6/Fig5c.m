close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% basic setting for fig5c
id.genotype = 'wt';
feedback_flags = [1, 0]; % no feedback and feedback
id.output ={'TNFmRNA','TNFpro','NFkBn','IKK','MyD88s','TRIFs','CpG', ...
            'CpGTLR9','PIC','PICTLR3','TLR4LPS'};
mRNA_filenames = {'./simData/mRNA_sim.csv','./simData/mRNA_sim_feedback.csv'};
sec_filenames = {'./simData/sec_sim.csv',['./simData/' ...
                    'sec_sim_feedback.csv']};
nfkbn_filenames = {'./simData/nfkb_sim_nofeedback.csv','./simData/nfkb_sim_feedback.csv'};

id.DT = 1;
id.sim_time = 480;
t = 0:id.sim_time;
plot_flag = 1 ;
stimuli.name= {'LPS','CpG','PIC','TNF'};

stimuli.dose= [100,100,50,1]; % 10ng/ml, 100nM, 50 ug/ml 

styles = {'k--','k-'};

for k = 1:2
    id.flag_noTnfFeedback = feedback_flags(k);%true or false
    
    for j = 1:3
        id.stimuli = stimuli.name{j}; %LPS, CpG, PIC,TNF. 
        id.dose = stimuli.dose(j); %'1','100' 
        
        sim{j} = getSimData(id);    
        % show the result 
        if plot_flag 
            
            % select figure
            if k ==1
                figure('position',[680   419   217*4   559/3*4])
            else
                figure(j)
            end
            
            for i=1:size(id.output,2)
                subplot(4,3,i);
                hold on
                plot(sim{j}(i,:),styles{k},'linewidth',2);
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
    tnfSec    = cumsum([sim{1}(2,:)*tp(9);sim{2}(2,:)*tp(9);sim{3}(2,:)*tp(9)/1.5]');
    nfkbn     = [sim{1}(3,:)/sim{1}(3,1);sim{2}(3,:)/sim{2}(3,1);sim{3}(3,:)/sim{3}(3,1)];

    csvwrite(mRNA_filenames{k},[t;tnfmRNA]')
    csvwrite(sec_filenames{k},[t;tnfSec']')
    csvwrite(nfkbn_filenames{k},[t;nfkbn]')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. plot in R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!R CMD BATCH Fig5c.R
!R CMD BATCH Fig6.R
!rm *.Rout
!rm *.RD* 
!rm *.Rhi*

