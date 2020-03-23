% Load in data and Speaker Index
load('K:\WorkingFolder\synchony\ot_synch_ff_del_fr_200_cor_spiketimes.mat')
load('K:\Data\Analysis\SpeakerIndex.mat','SpeakerIndex')

%% Run some quick calculations and look for the neuron pairs with strongest correlations

win = 10;

gm = sqrt(ot_synch.ff{3} .* ot_synch.ff{4});
synch = cellfun(@(x) ...
    sum( x(round(length(x)/2) - win : round(length(x)/2) + win)) / (win*2+1), ...
    ot_synch.ff{11});

ind_corr = diag(corr(gm', synch', 'rows', 'pairwise'));
     
[sorted, ix] = sort(ind_corr, 'descend');
ix(isnan(sorted)) = [];
ix = ix(1:10);

clear win

%% Plot just the top neuron pairs and look for good match between synch and fr

for p = 1:length(ix)
    figure(ix(p))
    subplot(1, 2, 1)
    SpRFplot_square(gm(ix(p), :)');
    title(num2str(ix(p)))
    subplot(1, 2, 2)
    SpRFplot_square(synch(ix(p), :)');    

end

clear p

%% Get the indices for particular az and el

good_pair = 133;
az = [-10, -30, -50];
el = 10;

ex = nan(1, 3);
for i = 1:length(az)
    idx = intersect(find(SpeakerIndex(:, 1) == az(i))', find(SpeakerIndex(:, 2) == el)');
    if ~isempty(idx)
        ex(i) = idx;
    else
        display(['no matching speaker for azimuth ' num2str(az(i)) ' and elevation ' num2str(el)])
    end
end

clear i idx az el 

%% Plot the full figure

figure
set(gcf, 'Position', [100 100 500 800])

subplot(4, 3, 1:6)
SpRFplot_square(gm(good_pair, :)');
title('Geometric firing rate spatial receptive field')

hold on

colors = ['b', 'y', 'r'];

for i = 1:length(ex)
    scatter(SpeakerIndex(ex(i), 1), SpeakerIndex(ex(i), 2), 100, colors(i), 'fill')
end

for i = 1:length(ex)
    subplot(4, 3, 6+i)
    n1 = ot_synch.ff{14}{good_pair}(ex(i), :)';
    n1 = cellfun(@(x) x*1000 - 50, n1, 'UniformOutput', 0);
    n2 = ot_synch.ff{15}{good_pair}(ex(i), :)';
    n2 = cellfun(@(x) x*1000 - 50, n2, 'UniformOutput', 0);
    y1 = get_spike_y(n1);
    y2 = get_spike_y(n2);
    n1 = n1(~cellfun('isempty',n1));
    n2 = n2(~cellfun('isempty',n2));
    scatter(cell2mat(n1), cell2mat(y1), 10, 'r', 'fill')
    hold on
    scatter(cell2mat(n2), cell2mat(y2), 10, 'b', 'fill')
    axis square
    xlabel('Time (ms)')
    ylabel('Repetition')
    ylim([0 20])
end

xx = -100:1:100;
xx = xx(2:end-1);

for i = 1:length(ex)
    subplot(4, 3, 9+i)
    plot(xx, ot_synch.ff{11}{good_pair, ex(i)}, colors(i))
    xlabel('Lag (ms)')
    ylabel('Normalized \newline coincidences')
    axis square
    ylim([0 0.001])
end

clear i n1 n2 y1 y2 xx colors