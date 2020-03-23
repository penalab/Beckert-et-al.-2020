%% Load in the data
load_repro_data

%% flatten the cell array

FR = cellfun(@(x) cell2mat(x'), FR, 'UniformOutput', 0);
REP = cellfun(@(x) cell2mat(x'), REP, 'UniformOutput', 0);

%% statistics

n_check = cellfun(@kstest, REP);

p_fr = ranksum(FR{1}, FR{2});
p_rep = ranksum(REP{1}, REP{2});

%%

boxplot2_points(1:2, REP(1:2), 0.75)
xlim([0 3])
xticks(1:2)
xticklabels({'itd', 'ild'})
ylabel('repro')
title('overall repro strength')

%%

bins = 0:5:100;

[Y_l, ~] = discretize(FR{2}, bins);
[Y_t, E] = discretize(FR{1}, bins);

itd = cell(length(E)-1, 1);
ild = cell(length(E)-1, 1);

for n = 1:length(E)-1
    itd{n} = REP{1}(Y_t == n);
    ild{n} = REP{2}(Y_l == n);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking one-to-one plots

load_repro_data

FR = cellfun(@(x) cell2mat(x'), FR, 'UniformOutput', 0);
REP = cellfun(@(x) cell2mat(x'), REP, 'UniformOutput', 0);

% FR = cellfun(@(x) round(x, 2), FR, 'UniformOutput', 0);


% Change these values to be what you want to compare
% dichotic would be 1 and 2
% free field would be 3 and 4
% ot ff would be 5 and 6
v1 = 2;
v2 = 1;

text = {'ITD', 'ILD', 'Az', 'El', 'otAz', 'otEl'};

remain = FR{v1};

rep_pairs = nan(length(FR{v2}), 2);
fr_pairs = nan(length(FR{v2}), 2);

for i = 1:length(FR{v2})
    
    difference = abs(remain - FR{v2}(i));
    
    rep_pairs(i, 2) = REP{v2}(i);
    fr_pairs(i, 2) = FR{v2}(i);
    
    [~, ix] = min(difference);
    
    rep_pairs(i, 1) = REP{v1}(ix);
    fr_pairs(i, 1) = FR{v1}(ix);
    
    % remain(ix) = [];
    
end

% % plotting

figure

subplot(2, 2, 1)
[r, p] = corr(fr_pairs(:, 1), fr_pairs(:, 2), 'rows', 'pairwise');
scatter_regression_plot(fr_pairs(:, 1), fr_pairs(:, 2))
xlabel([text{v1} ' firing rate'])
ylabel([text{v2} ' firing rate'])
title({'paired firing rate'; ['r = ' num2str(r) '; p = ' num2str(p)]})

subplot(2, 2, 2)
[r, p] = corr(rep_pairs(:, 1), rep_pairs(:, 2), 'rows', 'pairwise');
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