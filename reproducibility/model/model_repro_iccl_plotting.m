dir_name = 'K:\python\reproducibility\';

% itds = -200:20:200;
itds = 100;
SNR_ICcl = 2;

rates = cell(length(itds), 1);
repro = cell(length(itds), 1);

for t = 1:length(itds)
    load([dir_name '\data\iccl\repro_iccl_ITD' num2str(itds(t)) '_001_SNR' num2str(SNR_ICcl * 10) '.mat'], 'fr')
    load([dir_name '\data\iccl\repro_iccl_ITD' num2str(itds(t)) '_001_SNR' num2str(SNR_ICcl * 10) '_strength.mat'], 'strength')
    
    rates{t} = fr;
    repro{t} = strength;
    
end

rates = cell2mat(rates);
repro = cell2mat(repro);
    
clear dir_name t fr strength

%% quantification

rep_max = max(repro, [], 1);
rep_med = median(repro);

ind_corr = nan(size(rep_max));

for i = 1:size(ind_corr, 2)
    ind_corr(i) = corr(rates(:, i), repro(:, i), 'rows', 'pairwise');
end

clear i

%%

subplot(1, 3, 1)
boxplot2_points(1, {rep_max}, 0.75)
subplot(1, 3, 2)
boxplot2_points(1, {rep_med}, 0.75)
subplot(1, 3, 3)
boxplot2_points(1, {ind_corr}, 0.75)

%%

% f = norm01(rates, 2);
% r = norm01(repro, 2);

f = rates;
r = repro;

subplot(1, 10, 1:8)
% scatterheatmapplot(f(:), r(:), {'firing rate', 'reproducibility'}, 1)
scatter_regression_plot(f(:), r(:), 'r')
subplot(1, 10, 9:10)
boxplot2_points(1, {ind_corr}, 0.75)
xlim([0 2])
xticklabels('')