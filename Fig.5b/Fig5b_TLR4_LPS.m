
%% wt simulation
% 1. low dose, wt
id.dose = 10; %'1','100' 
id.DT = 1;
addpath('../src')
id.sim_time = 120
id.output ={'NFkBn','TNFnas', 'TNFmRNA','TNFpro','IKK'};
%'IKK','nfkb','irf'
genotypes ={'wt','mko','tko'};
plot_flag = 1;
for j = 1:3
    id.genotype = genotypes{j};
    sim{j} = getSimData(id);
    
    % show the result 
    if plot_flag 
        figure('position',[680   419   217*4   559/3*4])
        for i=1:size(id.output,2)
            subplot(4,3,i);
            plot(sim{j}(i,:),'k','linewidth',2);
            xlim([0 id.sim_time]);set(gca,'xtick',0:60:360);
            title(id.output{i})
        end
    end
end


%t(9) = .18;%tnf secretion rate , 2.5 fold less in tko

nfkbn     = [sim{1}(1,:);sim{2}(1,:);sim{3}(1,:)];
tnfNas    = [sim{1}(2,:);sim{2}(2,:);sim{3}(2,:)];
tnfmRNA   = [sim{1}(3,:);sim{2}(3,:);sim{3}(3,:)];
tnfPro    = [sim{1}(4,:);sim{2}(4,:);sim{3}(4,:)];
tnfSec    = cumsum([sim{1}(4,:)*.18;sim{2}(4,:)*.18;sim{3}(4,:)*.18/2.5]');
IKK    = [sim{1}(5,:);sim{2}(5,:);sim{3}(5,:)];
%%
t = 0:120;
csvwrite('./simData/nfkb_sim.csv',[t;nfkbn]')
csvwrite('./simData/nascent_sim.csv',[t;tnfNas]')
csvwrite('./simData/mRNA_sim.csv',[t;tnfmRNA]')
csvwrite('./simData/prot_sim.csv',[t;tnfPro]')
csvwrite('./simData/sec_sim.csv',[t;tnfSec']')

% end
!R CMD BATCH Fig5b_TLR4_LPS.R
