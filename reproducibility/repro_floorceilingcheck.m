load_repro_data

FR = cellfun(@(x) cell2mat(x'), FR, 'UniformOutput', 0);
REP = cellfun(@(x) cell2mat(x'), REP, 'UniformOutput', 0);

%%
   
text = {'iccl itd', 'iccl ild', 'iccl az', 'iccl el', 'ot az', 'ot el'};

for i = 1:6
    
subplot(3, 2, i)

fr = FR{i};
rep = REP{i};

[r, p] = corr(fr', rep', 'rows', 'pairwise');

scatter_regression_plot(fr, rep, 'r')
xlabel('Firing rate (Hz)')
ylabel('Reproducibility (coincidences/spike)')
title({['fr vs rep (' text{i} ')']; ['r = ' num2str(r) '; p = ' num2str(p)]})

end