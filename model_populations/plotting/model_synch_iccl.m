
ITD = 100;
rep = 1;
SNR = 2;
r_txt = sprintf('%03d', rep);

% dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\';
dir_name = 'K:\python\model_icx_codeonly\';

load([dir_name 'SynchResults\data\iccl\SynchICcl_ITD' num2str(ITD) '_' r_txt '_SNR' num2str(SNR * 10)], 'BITD', 'CF')
load([dir_name 'SynchResults\data\iccl_strength\SynchICcl_strength_ITD' num2str(ITD) '_' r_txt '_itd'])

%%

sim = nan(10000, 4);
dif = nan(10000, 4);

for f = 1:size(ccg_acrossitd, 2)
    for i = 1:size(ccg_acrossitd, 1)
        if abs(list_acrossitd(i, 1) - list_acrossitd(i, 2)) < 3
            L = find(isnan(sim(:, 1)), 1, 'first');
            sim(L, 1) = gm_acrossitd(i, f);
            sim(L, 2) = ccg_acrossitd(i, f);
            sim(L, 3) = sh_acrossitd(i, f);
            sim(L, 4) = cor_acrossitd(i, f);
        elseif abs(list_acrossitd(i, 1) - list_acrossitd(i, 2)) > 50
            L = find(isnan(dif(:, 1)), 1, 'first');
            dif(L, 1) = gm_acrossitd(i, f);
            dif(L, 2) = ccg_acrossitd(i, f);
            dif(L, 3) = sh_acrossitd(i, f);
            dif(L, 4) = cor_acrossitd(i, f);
        end
        
    end
end

L = find(isnan(sim(:, 1)), 1, 'first');
sim(L:end, :) = [];
L = find(isnan(dif(:, 1)), 1, 'first');
dif(L:end, :) = [];

%%

figure(1)
set(gcf, 'position', [0 0 1000 300])

subplot(1, 10, 1:7)
boxplot2_points(1:2, {sim(:, 2), dif(:, 2)}, 0.75)
boxplot2_points(4:5, {sim(:, 3), dif(:, 3)}, 0.75)
boxplot2_points(7:8, {sim(:, 4), dif(:, 4)}, 0.75)
xticks([1, 2, 4, 5, 7, 8])
xlim([0, 9])
ylabel('Synchrony (coincidences/spikes)')

%%

a = gm_acrossitd(:);
colors = {'k', 'r', 'b'};
subplot(1, 10, 8:10)
hold on

for t = 1:3
    switch t
        case 1
            b = ccg_acrossitd(:);
        case 2
            b = sh_acrossitd(:);
        case 3
            b = cor_acrossitd(:);
    end

scatter(a(1:100:end), b(1:100:end), 10, 'fill', colors{t})
xlabel('Geometric mean firing rate')
ylabel('Synchrony (coincidences/spike)')

% [cc, pp] = corr(a,b,'rows','pairwise');
% title(['r = ' num2str(cc) '; p = ' num2str(pp)])

end