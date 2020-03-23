dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\';
% dir_name = 'K:\python\model_icx_codeonly\';

load([dir_name 'SynchResults\model_out_full_mani.mat'])

%%

figure(2)
side_orig = nan(size(output_orig));
side_mani = nan(size(output_orig));

for t = 1:size(output_orig, 1)
    for r = 1:size(output_orig, 2)
        M = nanmax(output_orig{t, r}(bounds(t, 1):bounds(t, 2)));
        side_orig(t, r) = nanmax(output_orig{t, r}(bounds(t, 3):bounds(t, 4))) / M;
        side_mani(t, r) = nanmax(output_mani{t, r}(bounds(t, 3):bounds(t, 4))) / M;
    end
end

side_orig = side_orig(:)';
side_mani = side_mani(:)';

bar(1:2, [nanmean(side_orig, 2); nanmean(side_mani, 2)])

hold on

er = errorbar(1:2,[nanmean(side_orig); nanmean(side_mani)], [std(side_orig); std(side_mani)]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

xlim([0 3])
xticks(1:2)
xticklabels({'Initial', 'Manipulated'})
xtickangle(45)
ylabel({'Side-peak height', '% of main peak'})

hold off

%%
% peaks = bounds(:, 3) + round((bounds(:, 4) - bounds(:, 3)) / 2);
% 
% sidepeaks = nan(size(output_orig));
% 
% for t = 1:size(output_orig, 1)
%     for r = 1:size(output_orig, 2)
%         sidepeaks(t, r) = 1 - (output_orig{t, r}(peaks(t)) / output_mani{t, r}(peaks(t)));
%     end
% end
% 
% sidepeaks = sidepeaks(:);
% sidepeaks(isnan(sidepeaks)) = 0;
% 
% nanmean(sidepeaks)
% boxplot2_points(1, {sidepeaks}, 0.75);
% xlim([0 2])
% xticks(1)
% xticklabels('')
% ylabel('side-peak suppression')
%%

% sidepeaks = nan(size(output_orig));
% 
% for t = 1:size(output_orig, 1)
%     for r = 1:size(output_orig, 2)
%         o = output_orig{t, r}(bounds(t, 3):bounds(t, 4));
%         m = output_mani{t, r}(bounds(t, 3):bounds(t, 4));
%         sidepeaks(t, r) = 1 - (nanmax(o) / nanmax(m));
%     end
% end
% 
% sidepeaks = sidepeaks(:);
% sidepeaks(isnan(sidepeaks)) = 0;
% 
% nanmean(sidepeaks)
% boxplot2_points(1, {sidepeaks}, 0.75);
% xlim([0 2])
% xticks(1)
% xticklabels('')
% ylabel('side-peak suppression')
% 

%%
% peaks = bounds(:, 3) + round((bounds(:, 4) - bounds(:, 3)) / 2);
% 
% side_orig = nan(size(output_orig));
% side_mani = nan(size(output_orig));
% 
% for t = 1:size(output_orig, 1)
%     for r = 1:size(output_orig, 2)
%         side_orig(t, r) = nanmax(output_orig{t, r}(bounds(t, 3):bounds(t, 4)));
%         side_mani(t, r) = nanmax(output_mani{t, r}(bounds(t, 3):bounds(t, 4)));
%     end
% end
% 
% boxplot2_points(1:2, {side_orig(:), side_mani(:)}, 0.75);
% xlim([0 3])
% xticks(1:2)
% xticklabels({'original', 'manipulated'})
% ylabel('side-peak height')

%%
% itd = -200:20:200;
% 
% mani = cell(length(itd), 1);
% orig = cell(length(itd), 1);
% rate = cell(length(itd), 1);
% 
% for ITD = 1:length(itd)
%     
%     rep = 1;
%     r_txt = sprintf('%03d', rep);
% 
%     load([dir_name 'SynchResults\data\out_strength\SynchOut_strength_ITD' num2str(itd(ITD)) '_' r_txt])
% 
%     C = nan(1, nunits);
%     V = nan(1, nunits);
%     G = nan(1, nunits);
% 
%     C(1) = nanmean(syn_c{1}(1:2));
%     C(nunits) = nanmean(syn_c{1}(nunits-1:nunits));
%     V(1) = nanmean(syn_v{1}(1:2));
%     V(nunits) = nanmean(syn_v{1}(nunits-1:nunits));
%     G(1) = nanmean(gm(1:2));
%     G(nunits) = nanmean(gm(nunits-1:nunits));
% 
%     for i = 2:nunits-1
%         C(i) = nanmean(syn_c{1}(i, i-1:i+1));
%         V(i) = nanmean(syn_v{1}(i, i-1:i+1));
%         G(i) = nanmean(gm(i, i-1:i+1));
%     end
% 
%     mani{ITD} = C;
%     orig{ITD} = V;
%     rate{ITD} = G;
%     
% end
% 
% mani = cell2mat(mani);
% orig = cell2mat(orig);
% rate = cell2mat(rate);
% 
% %%
% synch_sidepeak = nan(size(orig, 1), 1);
% 
% for i = 1:length(peaks)
%     synch_sidepeak(i) = mani(i, peaks(i)) /orig(i, peaks(i));
% end
% 
% nanmean(synch_sidepeak)
% nanstd(synch_sidepeak)
% plot(synch_sidepeak)
% 
% %%
% peak_synch = nan(size(rate, 1), 1);
% 
% for i = 1:size(rate, 1)
%     [~, ix] = max(rate(i, :));
%     peak_synch(i) = orig(i, ix);
% end
% 
% %%
% 
% corr(orig(:), rate(:), 'rows', 'pairwise');