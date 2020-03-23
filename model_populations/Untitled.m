load('ICclFilter.mat')
%%

ccg_curve_orig = cellfun(@(x) x{5}{1}, synch_orig, 'UniformOutput', 0);
sh_curve_orig = cellfun(@(x) x{6}{1}, synch_orig, 'UniformOutput', 0);
cor_curve_orig = cellfun(@(x) x{7}{1}, synch_orig, 'UniformOutput', 0);

ccg_curve_filt = cellfun(@(x) x{5}{1}, synch_filt, 'UniformOutput', 0);
sh_curve_filt = cellfun(@(x) x{6}{1}, synch_filt, 'UniformOutput', 0);
cor_curve_filt = cellfun(@(x) x{7}{1}, synch_filt, 'UniformOutput', 0);

ccg_curve_orig = nanmean(cell2mat(ccg_curve_orig(:)'),2);
sh_curve_orig = nanmean(cell2mat(sh_curve_orig(:)'),2);
cor_curve_orig = nanmean(cell2mat(cor_curve_orig(:)'),2);

ccg_curve_filt = nanmean(cell2mat(ccg_curve_filt(:)'),2);
sh_curve_filt = nanmean(cell2mat(sh_curve_filt(:)'),2);
cor_curve_filt = nanmean(cell2mat(cor_curve_filt(:)'),2);
%%
figure
subplot(3,6,1:3)
hold on

scatter(gm(:), ccg_orig(:), 6, 'b')
scatter(gm(:), sh_orig(:), 3, 'r', 'fill')
scatter(gm(:), cor_orig(:), 3, 'k', 'fill')

lims = ylim;
ylim(lims)
xlabel('gm firing rate')
ylabel('synchrony')
title('original')

subplot(3,6,4:6)
hold on

scatter(gm(:), ccg_filt(:), 6, 'b')
scatter(gm(:), sh_filt(:), 3, 'r', 'fill')
scatter(gm(:), cor_filt(:), 3, 'k', 'fill')

ylim(lims)
xlabel('gm firing rate')
ylabel('synchrony')
title('filters')

N = size(ccg_orig, 1) * size(ccg_orig, 2);

subplot(3,6,7:8)
before_after_plot(ccg_orig(:), ccg_filt(:))
title(num2str(sum((ccg_orig(:) - ccg_filt(:)) > 0) / N))
ylabel('synchrony')
xlim([0 3])
xticks([1,2])
xticklabels({'orig', 'filters'})
subplot(3,6,9:10)
before_after_plot(sh_orig(:), sh_filt(:))
title(num2str(sum((sh_orig(:) - sh_filt(:)) > 0) / N))
xlim([0 3])
xticks([1,2])
xticklabels({'orig', 'filters'})
subplot(3,6,11:12)
before_after_plot(cor_orig(:), cor_filt(:))
title(num2str(sum((cor_orig(:) - cor_filt(:)) > 0) / N))
xlim([0 3])
xticks([1,2])
xticklabels({'orig', 'filters'})

X = linspace(-100, 100, 199)';
subplot(3,6,13:15)
hold on
plot(X, ccg_curve_orig, 'b', 'LineWidth', 3)
plot(X, sh_curve_orig, 'r')
plot(X, cor_curve_orig, 'k')
xlabel('time (ms)')
ylabel('coincidences')
title('original')
lims = ylim;
ylim(lims)

subplot(3,6,16:18)
hold on
plot(X, ccg_curve_filt, 'b', 'LineWidth', 3)
plot(X, sh_curve_filt, 'r')
plot(X, cor_curve_filt, 'k')
xlabel('time (ms)')
ylabel('coincidences')
title('filters')
ylim(lims)
