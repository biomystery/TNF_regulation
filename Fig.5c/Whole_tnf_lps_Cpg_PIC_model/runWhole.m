%% wt simulation
id.genotype = 'wt';
id.flag_noTnfFeedback = true;%true or false
% id.stimuli = 'TNF';
id.stimuli = 'CpG'; %LPS, CpG, PIC,TNF. 

% 1. low dose, wt
id.dose = 10; %'1','100' 
id.output ={'CpG','TLR4LPS','TLR4LPSen','MyD88s','TRIFs','IKKK','IKK',...
    'NFkBn','TNFt','TNF'};
sim{1} = getSimData(id);
%%
colors = {'k'};
figure('position',[680   419   217*4   559/3*4])
 
for i=1:size(id.output,2)
    subplot(4,3,i);
    plot(sim{1}(i,:),colors{mod(i,1)+1},'linewidth',2);
    xlim([0 360]);set(gca,'xtick',0:60:360);
    title(id.output{i})
end
