dir_name = 'K:\Code\Reproducibility\python\reproducibility\';

itds = -200:20:200;

rates = cell(length(itds), 1);
repro = cell(length(itds), 1);

for t = 1:length(itds)
    load([dir_name '\data\iccl\repro_iccl_ITD' num2str(itds(t)) '_001.mat'], 'fr')
    load([dir_name '\data\iccl\repro_iccl_ITD' num2str(itds(t)) '_001_strength.mat'], 'strength')
    
    rates{t} = fr;
    repro{t} = strength;
    
end

rates = cell2mat(rates);
repro = cell2mat(repro);
    
clear dir_name t fr strength

%% quantification

rep_max = max(repro, [], 1);

ind_corr = nan(size(rep_max));

for i = 1:size(ind_corr, 2)
    ind_corr(i) = corr(rates(:, i), repro(:, i), 'rows', 'pairwise');
end

clear i

%%

subplot(1, 2, 1)
boxplot2_points(1, {rep_max}, 0.75)
subplot(1, 2, 2)
boxplot2_points(1, {ind_corr}, 0.75)

%%

f = norm01(rates, 2);
r = norm01(repro, 2);

subplot(1, 10, 1:7)
scatterheatmapplot(f(:), r(:), {'Normalized firing rate', 'Normalized reproducibility'})
subplot(1, 10, 9:10)
boxplot2_points(1, {ind_corr}, 0.75)
xlim([0 2])
xticklabels('')
ylim([-1 1])
ylabel('Correlation coefficient (r)')