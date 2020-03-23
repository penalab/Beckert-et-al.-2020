
load('model_synch_icx_shifttest_20.mat')

n = 60;

subplot(2, 2, 1)
scatter(cell2mat(spikesICx(n, :)), cell2mat(get_spike_y(spikesICx(n, :))'), 5, 'r', 'fill');
hold on 
scatter(cell2mat(spikesICx(n+1, :)), cell2mat(get_spike_y(spikesICx(n+1, :))'), 5, 'b', 'fill');

xlabel('time (ms)')
ylabel('rep')
title('original')

for i = 1:size(spikesICx, 1)
    spikesICx(i, :) = spikesICx(i, randperm(size(spikesICx, 2)));
end

subplot(2, 2, 2)
scatter(cell2mat(spikesICx(n, :)), cell2mat(get_spike_y(spikesICx(n, :))'), 5, 'r', 'fill');
hold on 
scatter(cell2mat(spikesICx(n+1, :)), cell2mat(get_spike_y(spikesICx(n+1, :))'), 5, 'b', 'fill');

xlabel('time (ms)')
ylabel('rep')
title('shuffled')


load('model_synch_icx_shifttest_50.mat')

n = 60;

subplot(2, 2, 3)
scatter(cell2mat(spikesICx(n, :)), cell2mat(get_spike_y(spikesICx(n, :))'), 5, 'r', 'fill');
hold on 
scatter(cell2mat(spikesICx(n+1, :)), cell2mat(get_spike_y(spikesICx(n+1, :))'), 5, 'b', 'fill');

xlabel('time (ms)')
ylabel('rep')
title('original')

for i = 1:size(spikesICx, 1)
    spikesICx(i, :) = spikesICx(i, randperm(size(spikesICx, 2)));
end

subplot(2, 2, 4)
scatter(cell2mat(spikesICx(n, :)), cell2mat(get_spike_y(spikesICx(n, :))'), 5, 'r', 'fill');
hold on 
scatter(cell2mat(spikesICx(n+1, :)), cell2mat(get_spike_y(spikesICx(n+1, :))'), 5, 'b', 'fill');

xlabel('time (ms)')
ylabel('rep')
title('shuffled')