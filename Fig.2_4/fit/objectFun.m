function rmsd = objectFun(inputp)

residues = calScore(inputp); 
rmsd = sqrt(sum(residues.^2/numel(residues)));
end