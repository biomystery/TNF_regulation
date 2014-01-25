addpath('../src/')

sec_exp = csvread('../expdata/TNF_secrection.csv',1,0);
pro_exp = csvread('../expdata/combineTNF.csv',1,0);

fig4_proTNF
fig4_proTNF_sec
fig4_proTNF_tl
fig4_proTNF_both
!R CMD BATCH Fig4.R