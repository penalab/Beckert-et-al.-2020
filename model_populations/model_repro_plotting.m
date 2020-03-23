load('f:\desktop\Code\Reproducibility\python\reproducibility\data\out\repro_ITD100_002_constant_strength.mat')
c = strength;
load('f:\desktop\Code\Reproducibility\python\reproducibility\data\out\repro_ITD100_002_varying_strength.mat')
v = strength;
%%
load('f:\desktop\Code\Reproducibility\model_ICx\SynchResults\data\out\SynchOut_ITD100_002.mat')
sp_c = cellfun(@length, spikes2);
sp_v = cellfun(@length, spikesICx);
scatter(cell2mat(spikesICx(:)'), cell2mat(spikes2(:)'));

%%

corr(mean(sp_c, 2), mean(sp_v, 2))
corr(c', v', 'rows', 'pairwise')

figure(2)
subplot(2, 4, 5:6)
plot(bestITD, mean(sp_c, 2), 'LineWidth', 5)
hold on
plot(bestITD, mean(sp_v, 2), 'LineWidth', 2)
xlabel('best ITD')
ylabel('firing rate')
axis square
title('population response')
% subplot(3, 4, 9:10)
% scatter(c, v)
% xlabel('modulated ICx')
% ylabel('original ICx')
% title('reproducibility between populations')
subplot(2, 4, 7:8)
plot(bestITD, c)
hold on
plot(bestITD, v)
xlabel('best ITD')
ylabel('reproducibility')
title('reproducibility across population')
axis square

%%

spk = spikesICx;
spk2 = spikes2;

% ex = [15, 35, 65, 70];
ex = 1:size(spk, 1);

for i = 1:length(ex)
    
    %     subplot(2, 4, i)
    y = get_spike_y(spk(ex(i), :));
    scatter(cell2mat(spk(ex(i), :)), cell2mat(y), 10, 'fill')
    hold on
    y = get_spike_y(spk2(ex(i), :));
    scatter(cell2mat(spk2(ex(i), :)), cell2mat(y), 30)
    title(num2str(bestITD(ex(i))))
    axis([0 300 0 100])
    axis square
    pause
    close all
end