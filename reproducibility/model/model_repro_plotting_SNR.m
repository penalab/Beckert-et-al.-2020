<<<<<<< HEAD
load('K:\python\reproducibility\data\ICxSpikes_ITD0_SNR2.mat')
load('K:\python\reproducibility\data\repro_ITD0_SNR2.mat')
load('K:\python\reproducibility\data\repro_ITD0_SNR2_strength.mat')
=======
load('f:\desktop\Code\python\reproducibility\data\ICxSpikes_ITD0_tau.mat')
load('f:\desktop\Code\python\reproducibility\data\repro_ITD0_tau.mat')
load('f:\desktop\Code\python\reproducibility\data\repro_ITD0_tau_strength.mat')
>>>>>>> d36d3bb50474b38998a4bfa0967db02c3433a483

snr = reshape(repmat(tau, size(spikesICx, 1), 1), (size(tau, 2) * size(spikesICx, 1)), 1);

rep = cell(length(tau), 1);

for i = 1:length(tau)
    ix = snr == tau(i);
    rep{i} = [fr(ix)', strength(ix)'];
end

figure(1)
for i = 1:length(rep)
    subplot(2, 2, i)
    scatter(rep{i}(:, 1), rep{i}(:, 2), 10, 'fill')
    title(num2str(tau(i)))
    xlabel('firing rate')
    ylabel('reproducibility')
end

<<<<<<< HEAD
%%

load('K:\python\reproducibility\data\ICxSpikes_ITD0_tau.mat')
load('K:\python\reproducibility\data\repro_ITD0_tau.mat')
load('K:\python\reproducibility\data\repro_ITD0_tau_strength.mat')

t = reshape(repmat(tau, size(spikesICx, 1), 1), (size(tau, 2) * size(spikesICx, 1)), 1);

rep = cell(length(tau), 1);

for i = 1:length(tau)
    ix = t == tau(i);
    rep{i} = [fr(ix)', strength(ix)'];
end

for i = 1:length(rep)
    subplot(2, 2, i)
    scatter(rep{i}(:, 1), rep{i}(:, 2), 10, 'fill')
    title(num2str(tau(i)))
    xlabel('firing rate')
    ylabel('reproducibility')
=======


for c = 1:size(curve, 1)
    curve(c, :) = (curve(c, :) - fr(c) ^ 2) * cw / fr(c);
end

center = round(size(curve, 2) / 2);
win = (size(curve, 2) - center) * cw;
x = -win:cw:win;

figure(2)
plot(x, curve)

cur = cell(length(unique(params)), 1);
for p = unique(params)
    cur{p+1} = curve(params == p, :);
end

clear p win center

c = cellfun(@(x) nanmean(x, 1), cur, 'UniformOutput', 0);

figure(3)
for i = 1:length(cur)
    subplot(4, 1, i)
    plot(x, c{i})
    ylabel(num2str(tau(i)))
    xlim([-0.02, 0.02])
    axis square
>>>>>>> d36d3bb50474b38998a4bfa0967db02c3433a483
end