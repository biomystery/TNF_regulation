%% wt simulation 
% 1. low dose, wt
id.dose = 10; %'1','100' 
id.DT = 1;
%addpath('~/Dropbox/matlab/TLR4')

id.inputP  = [];%input.value;%[0,0,0,0];%input;  % input parameters
id.inputPid = [];%input.id;%[45:47,40];%[1:4, 9:18, 20:27, 35:48, 49:52];
id.inputvP = [];
id.inputvPid = [];

id.output ={'NFkBn','TNFnas', 'TNFmRNA','TNFpro','IKK'};%,'TNFsec'};
%'IKK','nfkb','irf'

id.genotype = 'wt';
sim{1} = getSimData(id);

id.genotype = 'mko';
sim{2} = getSimData(id);

id.genotype = 'tko';
sim{3} = getSimData(id);
%


%t(9) = .18;%tnf secretion rate 

nfkbn     = [sim{1}(1,:);sim{2}(1,:);sim{3}(1,:)];
tnfNas    = [sim{1}(2,:);sim{2}(2,:);sim{3}(2,:)];
tnfmRNA   = [sim{1}(3,:);sim{2}(3,:);sim{3}(3,:)];
tnfPro    = [sim{1}(4,:);sim{2}(4,:);sim{3}(4,:)];
tnfSec    = cumsum([sim{1}(4,:)*.18;sim{2}(4,:)*.18;sim{3}(4,:)*.18/1.5]);
IKK    = [sim{1}(5,:);sim{2}(5,:);sim{3}(5,:)];
%%
t = 0:240;
csvwrite('./simData/nfkb_sim.csv',[t;nfkbn]')
csvwrite('./simData/nascent_sim.csv',[t;tnfNas]')
csvwrite('./simData/mRNA_sim.csv',[t;tnfmRNA]')
csvwrite('./simData/prot_sim.csv',[t;tnfPro]')
csvwrite('./simData/sec_sim.csv',[t;tnfSec]')

% end
!R CMD BATCH Fig5b_TLR4_LPS.R
