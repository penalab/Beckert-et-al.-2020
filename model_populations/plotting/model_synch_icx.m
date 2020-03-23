itds = -200:20:200;
rep = 1;
r_txt = sprintf('%03d', rep);

dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\';

sim = nan(10000, 4);
dif = nan(10000, 4);

GM = cell(length(itds), 1);
SYN = cell(length(itds), 1);

for ITD = 1:length(itds)

load([dir_name 'SynchResults\data\out_strength\SynchOut_strength_ITD' num2str(itds(ITD)) '_' r_txt])

GM{ITD} = gm_var;
SYN{ITD} = syn_var;

[~, ix] = max(gm_var);
ix = synch_varying{ix}(1);

for i = 1:size(syn_var, 1)
    if synch_varying{i}(1) == ix || synch_varying{i}(2) == ix
    if abs(synch_varying{i}(1) - synch_varying{i}(2)) < 3
        L = find(isnan(sim(:, 1)), 1, 'first');
        sim(L, 1) = gm_var(i);
        sim(L, 2) = syn_var(i, 1);
        sim(L, 3) = syn_var(i, 2);
        sim(L, 4) = syn_var(i, 3);
    elseif abs(synch_varying{i}(1) - synch_varying{i}(2)) > 50
        L = find(isnan(dif(:, 1)), 1, 'first');
        dif(L, 1) = gm_var(i);
        dif(L, 2) = syn_var(i, 1);
        dif(L, 3) = syn_var(i, 2);
        dif(L, 4) = syn_var(i, 3);
    end
    end
end

end

L = find(isnan(sim(:, 1)), 1, 'first');
sim(L:end, :) = [];
L = find(isnan(dif(:, 1)), 1, 'first');
dif(L:end, :) = [];

%%

figure(1)
set(gcf, 'position', [0 0 1000 300])

subplot(1, 10, 1:7)
boxplot2_points(1:2, {sim(:, 2), dif(:, 2)}, 0.75)
boxplot2_points(4:5, {sim(:, 3), dif(:, 3)}, 0.75)
boxplot2_points(7:8, {sim(:, 4), dif(:, 4)}, 0.75)
xticks([1, 2, 4, 5, 7, 8])
xlim([0, 9])
ylabel({'Synchrony', '(coincidences/spikes)'})


%%

a = cell2mat(GM);
colors = {'k', 'r', 'b'};
SYN = cell2mat(SYN);

subplot(1, 10, 8:10)
hold on

for t = 1:3

    b = SYN(:, t);

    scatter(a(1:100:end), b(1:100:end), 10, 'fill', colors{t})
    xlabel('Geometric mean firing rate')
    ylabel({'Synchrony', '(coincidences/spike)'})

end