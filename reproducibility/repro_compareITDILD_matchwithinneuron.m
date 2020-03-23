%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking one-to-one plots

load_repro_data

% FR = cellfun(@(x) cell2mat(x'), FR, 'UniformOutput', 0);
% REP = cellfun(@(x) cell2mat(x'), REP, 'UniformOutput', 0);



% Change these values to be what you want to compare
% dichotic would be 1 and 2
% free field would be 3 and 4
% ot ff would be 5 and 6
v1 = 1;
v2 = v1+1;

text = {'ITD', 'ILD', 'Az', 'El', 'otAz', 'otEl'};

rep_pairs = nan(200, 2);
fr_pairs = nan(200, 2);

for n = 1:length(FR{v2})
    for i = 1:length(FR{v2}{n})
        
        difference = abs(FR{v1}{n} - FR{v2}{n}(i));
        
        L = find(isnan(rep_pairs(:, 1)), 1, 'first');
        
        rep_pairs(L, 2) = REP{v2}{n}(i);
        fr_pairs(L, 2) = FR{v2}{n}(i);
        
        [~, ix] = min(difference);
        
        rep_pairs(L, 1) = REP{v1}{n}(ix);
        fr_pairs(L, 1) = FR{v1}{n}(ix);
        
    end
end

L = find(isnan(rep_pairs(:, 1)), 1, 'first');
rep_pairs(L:end, :) = [];
fr_pairs(L:end, :) = [];

% % plotting

figure

subplot(2, 2, 1)
[r, p] = corr(fr_pairs(:, 1), fr_pairs(:, 2), 'rows', 'pairwise', 'type', 'Spearman');
scatter_regression_plot(fr_pairs(:, 1), fr_pairs(:, 2))
xlabel([text{v1} ' firing rate'])
ylabel([text{v2} ' firing rate'])
title({'paired firing rate'; ['r = ' num2str(r) '; p = ' num2str(p)]})

subplot(2, 2, 2)
[r, p] = corr(rep_pairs(:, 1), rep_pairs(:, 2), 'rows', 'pairwise', 'type', 'Spearman');
scatter_regression_plot(rep_pairs(:, 1), rep_pairs(:, 2))
xlabel([text{v1} ' reproducibility'])
ylabel([text{v2} ' reproducibility'])
title({[text{v1} ' vs ' text{v2} ' for matching firing rates']; ['r = ' num2str(r) '; p = ' num2str(p)]})

subplot(2, 2, 3)
before_after_plot(fr_pairs(:, 1), fr_pairs(:, 2))
xlim([0, 3])
xticks(1:2)
xticklabels({text{v1}, text{v2}})
ylabel('firing rate')



% % statistics

n_check = kstest(rep_pairs);

p_pairs = signrank(rep_pairs(:, 1), rep_pairs(:, 2));

mdl = fitlm(rep_pairs(:, 1), rep_pairs(:, 2), 'linear');
tbl = anova(mdl);

subplot(2, 2, 4)
before_after_plot(rep_pairs(:, 1), rep_pairs(:, 2))
xlim([0, 3])
xticks(1:2)
xticklabels({text{v1}, text{v2}})
ylabel('Reproducibility')
title(['p = ' num2str(p_pairs)])

%% 

d = rep_pairs(:, 1) - rep_pairs(:, 2);
h = kstest(d);

m = mean(d);
s = std(d);

p = signtest(d);


%%
fr_m = mean(fr_pairs, 2);
fr = fr_pairs(:, 1) - fr_pairs(:, 2);
D = rep_pairs(:, 1) - rep_pairs(:, 2);

subplot(2, 2, 1)
[r, p] = corr(fr, D, 'rows', 'pairwise');
scatter_regression_plot(fr, D)
xlabel('difference in firing rate')
ylabel('difference in reproducibility')
title({'Difference in Reproducibility (ITD v ILD)'; ['r = ' num2str(r) '; p = ' num2str(p)]})

subplot(2, 2, 2)
hold on
scatter_regression_plot(fr_m, rep_pairs(:, 1), 'r')
scatter_regression_plot(fr_m, rep_pairs(:, 2), 'b')
for i = 1:length(rep_pairs)
    plot([fr_m(i), fr_m(i)], rep_pairs(i, :), 'k')
end
xlabel('firing rate')
ylabel('reproducibility')

subplot(2, 2, 4)
boxplot2_points(1, {D}, 0.75)
xlim([0 2])
xticks([])
ylabel('difference in reproducibility')
title('ITD vs ILD')
axis square

subplot(2, 2, 3)
scatter_regression_plot(fr_m, D)
[r, p] = corr(fr_m, D, 'rows', 'pairwise');
xlabel('firing rate')
ylabel('difference in reproducibility')
title({'Difference in Reproducibility (ITD v ILD)'; ['r = ' num2str(r) '; p = ' num2str(p)]})

%%

fr = abs(fr_pairs(:, 1) - fr_pairs(:, 2));
rep = abs(rep_pairs(:, 1) - rep_pairs(:, 2));

scatter_regression_plot(fr, rep);
[r, p] = corr(fr, rep, 'rows', 'pairwise', 'type', 'Spearman');
xlabel('difference in firing rate')
ylabel('difference in reproducibility')
title({['r = ' num2str(r) '; p = ' num2str(p)]})
