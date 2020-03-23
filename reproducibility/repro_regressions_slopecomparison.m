%% Load in the data
load_repro_data

%% Calculate linear regressions for each neuron

n_terms = length(FR);

slopes = cell(n_terms, 1);
intercept = cell(n_terms, 1);

for n = 1:n_terms
    
    n_units = length(FR{n});
    slopes{n} = nan(n_units, 1);
    intercept{n} = nan(n_units, 1);
    
    for i = 1:n_units
        mdl = fitlm(FR{n}{i}', REP{n}{i}', 'linear');
        tbl = anova(mdl);
        if tbl.pValue(1) < 1 % Change this here if you only want significant variables
            slopes{n}(i) = mdl.Coefficients.Estimate(end);
            intercept{n}(i) = mdl.Coefficients.Estimate(1);
        end
    end
    
end

clear n_units n_terms n i mdl tbl

%% Quick plotting
% Here just relationship between firing rate of each neuron 
% (min left, max right) and the final slope

range_min = cell(size(FR));
range_max = cell(size(FR));
range_dyn = cell(size(FR));
for i = 1:length(FR)
    range_min{i} = cellfun(@nanmin, FR{i});
    range_max{i} = cellfun(@nanmax, FR{i});
    range_dyn{i} = range_max{i} - range_min{i};
end

figure(1)

for i = 1:length(range_min)
    subplot(6, 3, i*3-2)
    scatter_regression_plot(range_dyn{i}, slopes{i});
    subplot(6, 3, i*3-1)
    scatter_regression_plot(range_min{i}, slopes{i});
    subplot(6, 3, i*3)
    scatter_regression_plot(range_max{i}, slopes{i});
end

% And here the "before and after" plots between ITD/ILD or Az/El

figure(2)

for i = 1:3
    subplot(3, 3, i*3-2)
    before_after_plot(range_min{i*2-1}, range_min{i*2})
    xlim([0,3])
    xticks([1,2])
    xticklabels({'itd', 'ild'})
    subplot(3, 3, i*3-1)
    before_after_plot(range_max{i*2-1}, range_max{i*2})
    xlim([0,3])
    xticks([1,2])
    xticklabels({'itd', 'ild'})
    subplot(3, 3, i*3)
    before_after_plot(range_dyn{i*2-1}, range_dyn{i*2})
    xlim([0,3])
    xticks([1,2])
    xticklabels({'itd', 'ild'})
end

%% Statistics

% Check for normality
n_check_min = cellfun(@(x) kstest( (x-nanmean(x))/nanstd(x) ), range_min);
n_check_max = cellfun(@(x) kstest( (x-nanmean(x))/nanstd(x) ), range_max);
n_check_dyn = cellfun(@(x) kstest( (x-nanmean(x))/nanstd(x) ), range_dyn);
n_check_slopes = cellfun(@(x) kstest( (x-nanmean(x))/nanstd(x) ), slopes);

% Set up
p_nonpara_min = nan(3, 1);
p_para_min = nan(3, 1);
p_nonpara_max = nan(3, 1);
p_para_max = nan(3, 1);
p_nonpara_dyn = nan(3, 1);
p_para_dyn = nan(3, 1);
p_nonpara_slopes = nan(3, 1);
p_para_slopes = nan(3, 1);

% Calculate both parametric and non-parametric paired tests
% Confirm with the normality test which one you should use for reporting
for i = 1:3
    p_nonpara_min(i) = signrank(range_min{i*2-1}, range_min{i*2});
    [~, p_para_min(i)] = ttest(range_min{i*2-1}, range_min{i*2});
    p_nonpara_max(i) = signrank(range_max{i*2-1}, range_max{i*2});
    [~, p_para_max(i)] = ttest(range_max{i*2-1}, range_max{i*2});
    p_nonpara_dyn(i) = signrank(range_dyn{i*2-1}, range_dyn{i*2});
    [~, p_para_dyn(i)] = ttest(range_dyn{i*2-1}, range_dyn{i*2});
    p_nonpara_slopes(i) = signrank(slopes{i*2-1}, slopes{i*2});
    [~, p_para_slopes(i)] = ttest(slopes{i*2-1}, slopes{i*2});
end

% For reporting values if desired
reporting_min = [cellfun(@nanmean, range_min), cellfun(@nanstd, range_min)];
reporting_max = [cellfun(@nanmean, range_max), cellfun(@nanstd, range_max)];
reporting_dyn = [cellfun(@nanmean, range_dyn), cellfun(@nanstd, range_dyn)];
reporting_slopes = [cellfun(@nanmean, slopes), cellfun(@nanstd, slopes)];

clear i

%% Plot showing summary, just comparing ITD and ILD here

figure(3)
subplot(2, 2, 1)
boxplot2_points(1:2, slopes(1:2), 0.75)
xlim([0 3])
xticks(1:2)
xticklabels({'ITD', 'ILD'})
ylabel('slope')
title('slope comparison')

subplot(2, 2, 2)
boxplot2_points(1:2, range_min(1:2), 0.75)
xlim([0 3])
xticks(1:2)
xticklabels({'ITD', 'ILD'})
ylabel('FR')
title('min FR comparison')

subplot(2, 2, 3)
boxplot2_points(1:2, range_max(1:2), 0.75)
xlim([0 3])
xticks(1:2)
xticklabels({'ITD', 'ILD'})
ylabel('FR')
title('max FR comparison')

subplot(2, 2, 4)
boxplot2_points(1:2, range_dyn(1:2), 0.75)
xlim([0 3])
xticks(1:2)
xticklabels({'ITD', 'ILD'})
ylabel('FR')
title('dynamic range comparison')

%% Potential Figure Plots with Before and After Plots

figure(4)

subplot(2, 2, 1)
before_after_plot(slopes{1}, slopes{2})
xlim([0 3])
xticks(1:2)
xticklabels({'ITD', 'ILD'})
ylabel('slope')
title('slope comparison')

subplot(2, 2, 2)
before_after_plot(range_dyn{1}, range_dyn{2})
xlim([0 3])
xticks(1:2)
xticklabels({'ITD', 'ILD'})
ylabel('FR')
title('dynamic range comparison')

subplot(2, 2, 3)
before_after_plot(range_min{1}, range_min{2})
xlim([0 3])
xticks(1:2)
xticklabels({'ITD', 'ILD'})
ylabel('min firing rate')
title('min comparison')

subplot(2, 2, 4)
before_after_plot(range_min{1}, range_min{2})
xlim([0 3])
xticks(1:2)
xticklabels({'ITD', 'ILD'})
ylabel('FR')
title('max FR comparison')

