load('K:\WorkingFolder\synchony\ot_synch_ff_del.mat')

win = 10;

check = ot_synch.ff{6};
gm = sqrt(ot_synch.ff{3} .* ot_synch.ff{4});

% gm = norm01(gm);

for v = 11:13
    
    synch = cellfun(@(x) ...
        sum( x(round(length(x)/2) - win : round(length(x)/2) + win)) / (win*2+1), ...
        ot_synch.ff{v});
    synch(~check) = nan;
    
%     synch = norm01(synch);
    
    ind = diag(corr(gm', synch', 'rows', 'pairwise'));
    
    figure(v - 10)
    
    subplot(1, 10, 1:7)
    scatter_regression_plot(gm(:), synch(:))
    xlabel('normalized gm fr')
    ylabel('normalized synch')
    
    subplot(1, 10, 9:10)
    boxplot2_points(1, {ind}, 0.75)
    xticks(1)
    xticklabels([])
    ylabel('correlation coefficient')
    xlim([0 2])

end