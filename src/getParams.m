function pars = getParams()
% Use a hashtable or dictionary data structure. 
    params_name = {'V_tr','Km_tr','n','k_pr','kdeg_m','k_tl', ...
                   'kdeg_p','k_sec', 'pr_fold'};
    params_value = [11.2500,   0.5,     3, 0.6  , 0.02,   0.06,  ...
                   5.8e-2, 5.8e-2,    1.5];
    
    pars = containers.Map(params_name,params_value);
end