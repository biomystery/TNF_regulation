function pars = getParams()
% Use a hashtable or dictionary data structure. 
    params_name = {'V_tr','Km_tr','n','k_pr','kdeg_m','k_tl', ...
                   'kdeg_p','k_sec', 'pr_fold','tl_fold', 'sec_fold'};
    params_value = [1,   0.5,     3, 0.6  , 0.02,   0.06,  ...
                   0.06, 0.06,   1.5, 1.5, 5];
    
    pars = containers.Map(params_name,params_value);
end