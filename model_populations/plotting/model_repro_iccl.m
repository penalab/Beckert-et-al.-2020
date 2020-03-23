load('K:\python\model_icx_codeonly\SynchResults\data\iccl\SynchICcl_ITD160_001_SNR20.mat')
load('K:\python\reproducibility\data\iccl\repro_iccl_ITD160_001_SNR20.mat')

fr = repmat(fr, size(curve, 2), 1)';
curve = (curve - fr .^ 2) .* cw ./ fr;

%%

frs = 8;
itds = [50, 35, 25];

ticknum = 5;

figure
set(gcf, 'Position', [100 100 500 800])
subplot(4, 3, 1:6)

x = linspace(min(BITD), max(BITD), ticknum);
y = linspace(min(CF), max(CF), ticknum);

imagesc(mean(R, 3))

axis square
xticks(linspace(1, size(R, 2), ticknum))
xticklabels(x)
xlabel('Best ITD (µs)')
yticks(linspace(1, size(R, 1), ticknum))
yticklabels(y)
ylabel('Characteristic Frequency (Hz)')

colorbar

title('Model ICcl population response to 160 µs stimulus')

hold on

scatter(itds(1), frs, 100, 'fill', 'r')
scatter(itds(2), frs, 100, 'fill', 'y')
scatter(itds(3), frs, 100, 'fill', 'b')

subplot(4, 3, 7)
bsRaster(squeeze(spikesICcl(frs, itds(3), :)),[0 500])
ax = get(gca, 'Children');
set(ax, 'Color', 'b')
axis square
xlabel('Time (ms)')
ylabel('Repetition')
subplot(4, 3, 8)
bsRaster(squeeze(spikesICcl(frs, itds(2), :)),[0 500])
ax = get(gca, 'Children');
set(ax, 'Color', 'y')
axis square
xlabel('Time (ms)')
subplot(4, 3, 9)
bsRaster(squeeze(spikesICcl(frs, itds(1), :)),[0 500])
ax = get(gca, 'Children');
set(ax, 'Color', 'r')
axis square
xlabel('Time (ms)')

xx = -100:0.1:100;
xx = xx(2:end-1);

subplot(4, 3, 10)
rep_idx = find(itd == itds(3));
plot(xx, curve(rep_idx(frs), :), 'b')
ylabel('normalized coincidences')
xlabel('Lag (ms)')
ylim([-0.002, 0.015]);
axis square
subplot(4, 3, 11)
rep_idx = find(itd == itds(2));
plot(xx, curve(rep_idx(frs), :), 'y')
xlabel('Lag (ms)')
ylim([-0.002, 0.015]);
axis square
subplot(4, 3, 12)
rep_idx = find(itd == itds(1));
plot(xx, curve(rep_idx(frs), :), 'r')
xlabel('Lag (ms)')
ylim([-0.002, 0.015]);
axis square
%%


