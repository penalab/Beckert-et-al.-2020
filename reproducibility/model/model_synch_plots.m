group = 150;
fol = 'f:\desktop\WorkingFolder\model_brian';
load([fol '\ICxSpikes_ITD' num2str(group) '.mat'])
load([fol '\synch_model_ITD' num2str(group) '.mat'])
clear group fol
% %

win = 10;

gm = cellfun(@(x) sqrt(x{2} * x{3}), synch);
syn_ccg = cellfun(@(x) sum(x{4}{1}(round(length(x{4}{1})/2) - win : round(length(x{4}{1})/2) + win)) / (win*2+1), synch);
syn_sh = cellfun(@(x) sum(x{5}{1}(round(length(x{5}{1})/2) - win : round(length(x{5}{1})/2) + win)) / (win*2+1), synch);
syn_cor = cellfun(@(x) sum(x{6}{1}(round(length(x{6}{1})/2) - win : round(length(x{6}{1})/2) + win)) / (win*2+1), synch);

diff_itd = nan(length(pairs), 1);
diff_freq = nan(length(pairs), 1);

for p = 1:length(pairs)
    diff_itd(p) = abs(bestITD(pairs(p, 1)) - bestITD(pairs(p, 2)));
    diff_freq(p) = abs(bestF(pairs(p, 1)) - bestF(pairs(p, 2)));
end

tmp = syn_ccg;

describe = cell(7, 1);

describe{1} = max(tmp, [], 1);
describe{2} = nanmean(tmp, 1);
describe{3} = nan(1, size(gm, 2));
for p = 1:size(gm, 2)
    describe{3}(p) = nanmean(tmp(diff_itd <= 50, p), 1);
end
describe{4} = min(tmp, [], 1);
describe{5} = diag(corr(gm, tmp, 'rows', 'pairwise'), 0)';
describe{6} = corr(diff_itd, tmp, 'rows', 'pairwise');
describe{7} = corr(diff_freq, tmp, 'rows', 'pairwise');
describe = cell2mat(describe);

clear win p bestF bestITD pairs tmp
%%

col = distinguishable_colors(length(unique(gi)));
ue = unique(ge);

for e = 1:length(unique(ge))
    for i = 1:length(unique(gi))
        
    figure(e)
    subplot(3, 3, 1)
    hold on
    scatter(gm(:, e*5-5+i), syn_ccg(:, e*5-5+i), 10, col(i, :), 'fill');
    ylabel('standard ccg')
    subplot(3, 3, 2)
    hold on
    scatter(diff_itd, syn_ccg(:, e*5-5+i), 10, col(i, :), 'fill');
    title(['ge = ' num2str(ue(e))])
    subplot(3, 3, 3)
    hold on
    scatter(diff_freq, syn_ccg(:, e*5-5+i), 10, col(i, :), 'fill');
    subplot(3, 3, 4)
    hold on
    scatter(gm(:, e*5-5+i), syn_ccg(:, e*5-5+i), 10, col(i, :), 'fill');
    ylabel('shifted ccg')
    subplot(3, 3, 5)
    hold on
    scatter(diff_itd, syn_sh(:, e*5-5+i), 10, col(i, :), 'fill');
    subplot(3, 3, 6)
    hold on
    scatter(diff_freq, syn_sh(:, e*5-5+i), 10, col(i, :), 'fill');
    subplot(3, 3, 7)
    hold on
    scatter(gm(:, e*5-5+i), syn_cor(:, e*5-5+i), 10, col(i, :), 'fill');
    ylabel('corrected ccg')
    xlabel('gm firing rate')
    subplot(3, 3, 8)
    hold on
    scatter(diff_itd, syn_cor(:, e*5-5+i), 10, col(i, :), 'fill');
    xlabel('diff ITD')
    subplot(3, 3, 9)
    hold on
    scatter(diff_freq, syn_cor(:, e*5-5+i), 10, col(i, :), 'fill');
    xlabel('diff freq')
    legend(num2str(unique(gi)))
    end
end
clear col e i

%%

figure(4)
subplot(1, 3, 1)
imagesc(reshape(describe(2,:), 5, 3))
xticks(1:3)
xticklabels(unique(ge))
xlabel('ge')
yticks(1:5)
yticklabels(unique(gi))
ylabel('gi')
h = colorbar;
ylabel(h, 'synchrony')
title('mean synchrony')
axis square

subplot(1, 3, 3)
imagesc(reshape(describe(1,:), 5, 3))
xticks(1:3)
xticklabels(unique(ge))
xlabel('ge')
yticks(1:5)
yticklabels(unique(gi))
ylabel('gi')
h = colorbar;
ylabel(h, 'synchrony')
title('max synchrony')
axis square

subplot(1, 3, 2)
imagesc(reshape(describe(3,:), 5, 3))
xticks(1:3)
xticklabels(unique(ge))
xlabel('ge')
yticks(1:5)
yticklabels(unique(gi))
ylabel('gi')
h = colorbar;
ylabel(h, 'synchrony')
title('mean synchrony for diff ITD <= 50µs')
axis square

clear h
%%

for params = 1:length(ge)
    X = spikesICx(:, :, params);
    Y = get_spike_y(X);
    for n = 1:size(X, 1)
        x = cell2mat(X(n, :));
        y = cell2mat(Y(n, :));
        if ~isempty(x)
            figure(10)
            scatter(x, y, 10, 'fill')
            axis([0 300 0 100])
            xlabel('time (ms)')
            ylabel('reps')
            title([num2str(params) ' --- ' num2str(n)])
            pause
            close
        end
    end
end

clear params X Y n x y
%%
load('f:\desktop\WorkingFolder\synchony\ot_synch_ff_del_fr_200_cor.mat')

win = 10;
ot = cell(4, 1);

for i = 1:3
    ot{i+1} = cellfun(@(x) ...
        sum( x(round(length(x)/2) - win : round(length(x)/2) + win)) / (win*2+1), ...
        ot_synch.ff{10 + i});
end

ot{1} = sqrt(ot_synch.ff{3} .* ot_synch.ff{4});

ot_syn = nan(size(ot{1}, 1), 4);
[ot_syn(:, 1), idx] = max(ot{1}, [], 2);

for i = 1:size(ot{1})
    ot_syn(i, 2) = ot{2}(i, idx(i));
    ot_syn(i, 3) = ot{3}(i, idx(i));
    ot_syn(i, 4) = ot{4}(i, idx(i));
end

figure(5)
for i = 2:4
    subplot(1, 3, i-1)
    scatter(ot_syn(:, 1), ot_syn(:, i), 10, 'fill')
    if i == 2
        ylabel('standard ccg')
    elseif i == 3
        ylabel('shifted ccg')
        title('OT synchrony')
    else
        ylabel('corrected ccg')
    end
    xlabel('gm firing rate')
    ylim([-0.0005, 0.0035])
end

describe_ot = nanmean(ot_syn);

clear ot_synch win i idx ot