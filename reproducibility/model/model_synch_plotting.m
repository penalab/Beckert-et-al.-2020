load('f:\desktop\Code\Reproducibility\python\reproducibility\data\icxsynch_ITD150_2_spike1.mat')

gm = sqrt(icx_synch{2} .* icx_synch{3});

win = 10;

synch = cellfun(@(x) ...
    sum( x(round(length(x)/2) - win : round(length(x)/2) + win)) / (win*2+1), ...
    icx_synch{7});

subplot(1, 2, 1)
scatter(gm, synch, 50, 'fill')
title(num2str(corr(gm', synch', 'rows', 'pairwise')))
xlabel('gm firing rate')
ylabel('standard ccg')

load('f:\desktop\WorkingFolder\synchony\ot_synch_ff_del_fr.mat')

gm = sqrt(ot_synch.ff{3} .* ot_synch.ff{4});

win = 10;

synch = cellfun(@(x) ...
    sum( x(round(length(x)/2) - win : round(length(x)/2) + win)) / (win*2+1), ...
    ot_synch.ff{13});

subplot(1, 2, 2)
scatter(gm(:), synch(:), 10, 'fill', 'r')
title(num2str(corr(gm(:), synch(:), 'rows', 'pairwise')))
xlabel('gm firing rate')
ylabel('standard ccg')