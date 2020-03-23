ITD = 80;
rep = 1;
r_txt = sprintf('%03d', rep);

dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\';
% dir_name = 'K:\python\model_icx_codeonly\';

load([dir_name 'SynchResults\data\out_mani\SynchOut_strength_ITD' num2str(ITD) '_' r_txt])

load([dir_name 'SynchResults\data\iccl\SynchICcl_ITD' num2str(ITD) '_' r_txt])

load([dir_name 'SynchResults\model_out_full_mani.mat'])
output_mani = output_mani(find(itd == ITD, 1, 'first'), :);
output_orig = output_orig(find(itd == ITD, 1, 'first'), :);

output_mani = [mean(cell2mat(output_mani')); std(cell2mat(output_mani'))];
output_orig = [mean(cell2mat(output_orig')); std(cell2mat(output_orig'))];

clear ITD rep r_txt dir_name

% %
% close all

trial = 50;

figure(1)
set(gcf, 'Position', [0, 0, 800, 1200])
subplot(4, 6, 7:9)
% y = cell2mat(get_spike_y(spikesICx(:, trial))');
x = spikesICx(:, trial)';
y = cell2mat(cellfun(@(i, x) ones(length(x), 1) * i, num2cell(1:length(x)), x, 'UniformOutput', 0)');
y = bestITD(y);
scatter(cell2mat(spikesICx(:, trial)'), y, 5, 'b','fill')
xlabel('Time (ms)','fontsize',15)
ylabel('ITD (탎)','fontsize',15)
title('Initial','fontsize',15)

subplot(4, 6, 10:12)
% y = cell2mat(get_spike_y(spikes2(:, trial))');
x = spikes2(:, trial)';
y = cell2mat(cellfun(@(i, x) ones(length(x), 1) * i, num2cell(1:length(x)), x, 'UniformOutput', 0)');
y = bestITD(y);
scatter(cell2mat(spikes2(:, trial)'), y, 5, 'r', 'fill')
xlabel('Time (ms)','fontsize',15)
yticklabels('')
title('Manipulated','fontsize',15)



for t = 1:3

names = {'Standard', 'Shifted', 'Corrected'};

C = nanmean([[nan; diag(syn_c{t}, -1)],[diag(syn_c{t}, 1); nan]], 2);
V = nanmean([[nan; diag(syn_v{t}, -1)],[diag(syn_v{t}, 1); nan]], 2);
G = nanmean([[nan; diag(gm, -1)],[diag(gm, 1); nan]], 2);

% subplot(4, 2, t*2+3:t*2+4)
subplot(4, 6, t*2+11:t*2+12)
plot(bestITD, V, 'b', 'LineWidth', 3)
hold on
plot(bestITD, C, 'r', 'LineWidth', 1)
xlabel('ITD (탎)','fontsize',15)
if t == 1
    ylabel({'Synchrony', '(coincidences/spike)'},'fontsize',15)
else
    ylabel('')
    yticklabels('')
end
title([names{t} ' CCG'],'fontsize',15)
ylim([0 0.015])

end

subplot(4, 6, 2:5)
plot(bestITD, G, 'k', 'LineWidth', 3)
xlabel('ITD (탎)')
ylabel({'Geometric mean', 'firing rate'},'fontsize',15)
title('Space-map population response','fontsize',15)

subplot(4, 6, 20:23)
errorbar(bestITD, output_orig(1, :), output_orig(2, :), 'b', 'LineWidth', 3)
hold on
errorbar(bestITD, output_mani(1, :), output_mani(2, :), 'r', 'LineWidth', 1)
xlabel('ITD (탎)','fontsize',15)
ylabel('Spike count','fontsize',15)
axis square
title('Readout population response','fontsize',15)
ylim([-1 13])
% plot(bestITD, mean(squeeze(SC(1, :, :))))
% hold on
% plot(bestITD, mean(squeeze(SC2(1, :, :))))
%%

figure(2)

subplot(3, 1, 1)
imagesc(BITD,CF,Rm);axis xy;colorbar
xlabel('Best ITD (\mus)','fontsize',15)
ylabel('Best Frequency (Hz)','fontsize',15)
title('ICcl Spike Counts','fontsize',15)

subplot(3, 1, 2)
plot(bestITD, G, 'k', 'LineWidth', 3)
xlabel('ITD (탎)','fontsize',15)
ylabel('Geometric mean firing rate','fontsize',15)
title('Tuning curve')

subplot(3, 1, 3)
errorbar(bestITD, output_orig(1, :), output_orig(2, :), 'b', 'LineWidth', 3)
hold on
errorbar(bestITD, output_mani(1, :), output_mani(2, :), 'r', 'LineWidth', 1)
xlabel('ITD (탎)','fontsize',15)
ylabel('Spike count','fontsize',15)
title('Output population')

