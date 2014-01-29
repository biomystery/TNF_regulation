pars = getParams();

k  = pars.keys;
v = pars.values;
n = numel(k);

fid = fopen('pars.csv','w');
for i = 1: n
    fprintf(fid,'%d,%s,%4.2f\n',i,k{i},v{i});
end
fclose(fid);