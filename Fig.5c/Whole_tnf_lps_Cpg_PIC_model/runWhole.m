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