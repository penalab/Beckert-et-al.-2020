% current
load('K:\WorkingFolder\synchony\ot_synch_ff_del_fr_200_cor.mat')

win = 10;
norm = 0;
% bins = 0:0.1:1;
bins = 100;

for v = 11:13
    
    check = ot_synch.ff{6};
    gm = sqrt(ot_synch.ff{3} .* ot_synch.ff{4});
    % gm(~check) = nan;
    if norm
        gm = norm01(gm);
    end
    
    synch = ot_synch.ff{v};
    %     synch(~check) = {nan(199, 1)};
    
    gm_bin = nan(size(gm, 1), bins);
    syn_bin = cell(size(synch, 1), bins);
    
    for r = 1:size(gm_bin, 1)
        [~, ~, idx] = histcounts(gm(r, :), bins);
        for c = 1:size(gm_bin, 2)
            gm_bin(r, c) = nanmean(gm(r, idx == c));
            syn_bin{r, c} = nansum(cell2mat(synch(r, idx == c)), 2) / sum(idx == c);
            if isempty(syn_bin{r, c})
                syn_bin{r, c} = nan(199, 1);
            end
        end
    end
    
    synch = cellfun(@(x) ...
        sum( x(round(length(x)/2) - win : round(length(x)/2) + win)) / (win*2+1), ...
        syn_bin);
%     synch = cellfun(@max, syn_bin);
    
    if norm
        synch = norm01(synch);
    end
    
    gm = gm_bin;
    
    ind = diag(corr(gm', synch', 'rows', 'pairwise'));
    
    figure(100 + (v - 10))
    
    subplot(1, 10, 1:7)
    hold on
       
    a = gm(:);
    b = synch(:);
    
    a = a(isfinite(b));
    b = b(isfinite(b));
    
    [sl, inter] = scatter_regression_plot(a, b, 'r');
    disp(num2str(sl))
    disp(inter)
    
    % We can still add the axis labels in either condition
    xlabel('Geometric Mean Firing Rate')
    ylabel('Synchrony')
    
    % Run a quick Pearson's correlation to see the relationship and
    % significance, put it as the title for easy visualization
    [cc, pp] = corr(a,b,'rows','pairwise');
    title(['r = ' num2str(cc) '; p = ' num2str(pp)])
    
    subplot(1, 10, 9:10)
    boxplot2_points(1, {ind}, 0.75)
    xticks(1)
    xticklabels([])
    ylabel('Correlation coefficient (r)')
    xlim([0 2])
    ylim([-1 1])
    
    % and some stats
    
    Z = Fisher_Z(ind);
    disp(num2str(signrank(Z)))
    m = nanmean(ind);
    s = nanstd(ind);
    disp(['mean = ' num2str(m) '; std = ' num2str(s)])
end