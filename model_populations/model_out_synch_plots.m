ITD = 60;
rep = 1;
r_txt = sprintf('%03d', rep);

dir_name = 'K:\python\model_icx_codeonly\';

load([dir_name 'SynchResults\data\out_strength\SynchOut_strength_ITD' num2str(ITD) '_' r_txt])

clear ITD rep r_txt dir_name

%% Plot synchrony of individual neurons
t = 1;

for i = 1:nunits
    subplot(1, 2, 1)
    plot(syn_c{t}(i, :))
    title(num2str(i))
    subplot(1, 2, 2)
    plot(syn_v{t}(i, :))
    title(num2str(i))
    pause
end

%% Plot scatter of synchrony vs firing rate

%figure(1)
subplot(3, 2, 2)
scatter(gm_con, syn_con(:, 1), 10, 'r', 'fill')
title(['r = ' num2str(corr(gm_con, syn_con(:, 1), 'rows', 'pairwise'))])
xlabel('gm firing rate')
ylabel('standard synchrony')
axis square
subplot(3, 2, 1)
scatter(gm_var, syn_var(:, 1), 10, 'b', 'fill')
title(['r = ' num2str(corr(gm_var, syn_var(:, 1), 'rows', 'pairwise'))])
xlabel('gm firing rate')
ylabel('standard synchrony')
axis square
   
%% Calculate the synchrony of nearby neurons and plot result

C = nan(1, nunits);
V = nan(1, nunits);
G = nan(1, nunits);

C(1) = nanmean(syn_c{1}(1:2));
C(nunits) = nanmean(syn_c{1}(nunits-1:nunits));
V(1) = nanmean(syn_v{1}(1:2));
V(nunits) = nanmean(syn_v{1}(nunits-1:nunits));
G(1) = nanmean(gm(1:2));
G(nunits) = nanmean(gm(nunits-1:nunits));

for i = 2:nunits-1
    C(i) = nanmean(syn_c{1}(i, i-1:i+1));
    V(i) = nanmean(syn_v{1}(i, i-1:i+1));
    G(i) = nanmean(gm(i, i-1:i+1));
end

% figure(2)
subplot(3, 2, 3:4)
plot(bestITD, V, 'b', 'LineWidth', 3)
hold on
plot(bestITD, C, 'r', 'LineWidth', 1)
xlabel('best ITD')
ylabel('standard synchrony')
title('comparing manipulated data - synchrony')

%% Plot relationship of firing rate and synchrony of each neuron and their nearby neurons
limits = [-50, 50];

for i = 1:length(limits)
    idx = find(bestITD == limits(i));
    limits(i) = idx;
end

% figure(3)
subplot(3, 2, 5)
scatter(G(1:limits(1)), V(1:limits(1)), 10, 'b', 'fill')
hold on
scatter(G(1:limits(1)), C(1:limits(1)), 10, 'r', 'fill')
xlabel('gm firing rate')
ylabel('standard synchrony')
title('relationship of synchrony in side peak')
axis square
subplot(3, 2, 6)
scatter(G(limits(2):end), V(limits(2):end), 10, 'k', 'fill')
xlabel('gm firing rate')
ylabel('standard synchrony')
title('relationship of synchrony in main peak')
axis square


