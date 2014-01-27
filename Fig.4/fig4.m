addpath('../src/')

sec_exp = csvread('../expdata/TNF_secrection.csv',1,0);
pro_exp = csvread('../expdata/combineTNF.csv',1,0);

fig4_proTNF % no regulation 
fig4_proTNF_sec % sec regulation 
fig4_proTNF_tl % tl regulation 
fig4_proTNF_both % both regulation 
!R CMD BATCH Fig4.R