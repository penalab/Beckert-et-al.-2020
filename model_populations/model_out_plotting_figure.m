cd('K:\python\model_icx_codeonly\SynchResults\data\out_synch_multireps')

files = ls;
files = files(3:end, :);

syn_curve_orig_ccg = nan(size(files, 1), 81);
syn_curve_mani_ccg = nan(size(files, 1), 81);
syn_curve_orig_sh = nan(size(files, 1), 81);
syn_curve_mani_sh = nan(size(files, 1), 81);

spiketimes = cell(size(files, 1), 2);

out_resp = cell(size(files, 1), 2);

sidepeaks = nan(size(files, 1), 2);

for f = 1:size(files, 1)
    
    load(files(f, :), 'SC', 'SC2', 'syn_v', 'syn_c', 'ind2', 'spikesICx', 'spikes2')
    
    spiketimes{f, 1} = spikesICx;
    spiketimes{f, 2} = spikes2;
    
    M = mean(SC,1);
    M2 = mean(SC2,1);
    sidepeaks(f, 1) = max(M(ind2))/max(M);
    sidepeaks(f, 2) = max(M2(ind2))/max(M2);

    
    out_resp{f, 1} = M;
    out_resp{f, 2} = M2;
    
    syn_curve_orig_ccg(f, 1) = nanmean(syn_v{1}(1, 1:2));
    syn_curve_mani_ccg(f, 1) = nanmean(syn_c{1}(1, 1:2));
    syn_curve_orig_ccg(f, end) = nanmean(syn_v{1}(end, end-1:end));
    syn_curve_mani_ccg(f, end) = nanmean(syn_c{1}(end, end-1:end));
    
    for p = 2:80
        syn_curve_orig_ccg(f, p) = nanmean(syn_v{1}(p, p-1:p+1));
        syn_curve_mani_ccg(f, p) = nanmean(syn_c{1}(p, p-1:p+1));
    end
    
    syn_curve_orig_sh(f, 1) = nanmean(syn_v{2}(1, 1:2));
    syn_curve_mani_sh(f, 1) = nanmean(syn_c{2}(1, 1:2));
    syn_curve_orig_sh(f, end) = nanmean(syn_v{2}(end, end-1:end));
    syn_curve_mani_sh(f, end) = nanmean(syn_c{2}(end, end-1:end));
    
    for p = 2:80
        syn_curve_orig_sh(f, p) = nanmean(syn_v{2}(p, p-1:p+1));
        syn_curve_mani_sh(f, p) = nanmean(syn_c{2}(p, p-1:p+1));
    end  
    
end

spiketimes{1,1} = cellfun(@(x) x - 150, spiketimes{1,1}, 'UniformOutput', 0);
spiketimes{1,2} = cellfun(@(x) x - 150, spiketimes{1,2}, 'UniformOutput', 0);

out_resp_orig = cell2mat(out_resp(:, 1));
out_resp_mani = cell2mat(out_resp(:, 2));

itds = linspace(-200, 200, 81);
popresp = cellfun(@length, spiketimes{1, 1});

sidepeaks_synch(:, 1) = max(syn_curve_orig_ccg(:, ind2), [], 2);
sidepeaks_synch(:, 2) = max(syn_curve_mani_ccg(:, ind2), [], 2);
sidepeaks_synch(:, 3) = max(syn_curve_orig_sh(:, ind2), [], 2);
sidepeaks_synch(:, 4) = max(syn_curve_mani_sh(:, ind2), [], 2);

[~, p_synch_ccg] = ttest(sidepeaks_synch(:, 1), sidepeaks_synch(:, 2));
[~, p_synch_sh] = ttest(sidepeaks_synch(:, 3), sidepeaks_synch(:, 4));

sidepeaks_synch_report = [mean(sidepeaks_synch); std(sidepeaks_synch)];

sidepeaks_synch_perc = sidepeaks_synch_report(1,2) / sidepeaks_synch_report(1,1);

clear M M2 f p SC SC2 syn_v syn_c ind2 spikesICx spikes2 files out_resp


%%

subplot(2, 2, 1:2)
errorbar(itds, mean(popresp, 2), std(popresp, [], 2), 'k', 'LineWidth', 2)
xlabel('ITD (탎)')
ylabel('Spike Count')

% subplot(4, 2, 3);
% bsRaster(spiketimes{1, 1}(:, 50),[0 100]);hold on
% set(gca,'YTick',1:10:81,'YTickLabel',itds(1:10:81))
% xlabel('Time (ms)','fontsize',15)
% ylabel('Best ITD (\mus)','fontsize',15)
% title('Initial ICx spikes','fontsize',15)

% subplot(4, 2, 4);
% bsRaster(spiketimes{1, 2}(:, 50),[0 100]);hold on
% set(gca,'YTick',1:10:81,'YTickLabel',itds(1:10:81))
% xlabel('Time (ms)','fontsize',15)
% ylabel('Best ITD (\mus)','fontsize',15)
% title('Manipulated ICx spikes','fontsize',15)

% subplot(4, 2, 5)
% errorbar(itds, mean(syn_curve_orig_ccg, 1), std(syn_curve_orig_ccg, [], 1), 'b', 'LineWidth', 2)
% hold on
% errorbar(itds, mean(syn_curve_mani_ccg, 1), std(syn_curve_mani_ccg, [], 1), 'r')
% xlabel('ITD (탎)')
% ylabel('Synchrony (coindences/spike)')
% title('Standard CCG')

% subplot(4, 2, 6)
% errorbar(itds, mean(syn_curve_orig_sh, 1), std(syn_curve_orig_sh, [], 1), 'b', 'LineWidth', 2)
% hold on
% errorbar(itds, mean(syn_curve_mani_sh, 1), std(syn_curve_mani_sh, [], 1), 'r')
% xlabel('ITD (탎)')
% ylabel('Synchrony (coindences/spike)')
% title('Shifted CCG')

subplot(2, 2, 3)
errorbar(itds, mean(out_resp_orig, 1), std(out_resp_orig, [], 1), 'b', 'LineWidth', 2)
hold on
errorbar(itds, mean(out_resp_mani, 1), std(out_resp_mani, [], 1), 'r')
xlabel('ITD (탎)')
ylabel('Spike count')
title('Readout population response')

subplot(2, 2, 4)
bar(1:2, mean(sidepeaks))                
hold on
er = errorbar(1:2,mean(sidepeaks),std(sidepeaks));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
xticks(1:2)
xlim([0, 3])
xticklabels({'Initial', 'Manipulated'})
ylabel('Side peak height (% of main peak)')

[~, p_sidepeaks] = ttest(sidepeaks(:, 1), sidepeaks(:, 2));
title(num2str(p_sidepeaks))

clear er

sidepeaks_report = [mean(sidepeaks); std(sidepeaks)];
