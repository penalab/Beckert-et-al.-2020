load('model_synch_iccl_shifttest_100_shortdur.mat')
%%

neu = 1;

figure
subplot(3,2,1)
hold on

for i = 1:Nf
    plot(synch_acrossitd_orig{neu, i}{5}{1})
end
title('original')

subplot(3,2,2)
hold on

for i = 1:Nf
    plot(synch_acrossitd_shift{neu, i}{5}{1})
end
title('shifted')

%%

win = 10;

center = round(length(synch_acrossitd_shift{neu, 2}{5}{1}) / 2);
s = sum(synch_acrossitd_shift{neu, 2}{5}{1}(center - win:center + win)) / (win*2+1)
o = sum(synch_acrossitd_orig{neu, 2}{5}{1}(center - win:center + win)) / (win*2+1)

%%

subplot(3, 2, 3)
scatter(ccg_acrossitd_orig(:),ccg_acrossitd_shift(:), 5, 'r')
hold on
plot([0, 0.015], [0, 0.015], 'b')
title('original vs. shifted')
xlabel('original')
ylabel('shifted')

%%
%%
load('K:\WorkingFolder\synchony\ot_synch_ff_del_fr_200_cor.mat')

win = 10;
norm = 0;
% bins = 0:0.1:1;
bins = 100;

v = 12;

gm = sqrt(ot_synch.ff{3} .* ot_synch.ff{4});
synch = cellfun(@(x) ...
        sum( x(round(length(x)/2) - win : round(length(x)/2) + win)) / (win*2+1), ...
        ot_synch.ff{v});
    
subplot(3,2,6)
scatter_regression_plot(gm(:), synch(:))
xlabel('geometric mean')
ylabel('synch')
title('Physiological OT')
ylim([0, 0.0025])
xlim([0, 60])

%%
t = 2;

a = cell2mat(GM);
colors = {'k', 'r', 'b'};
syn = cell2mat(SYN);

subplot(3,2,5)

b = syn(:, t);

scatter_regression_plot(a(1:end), b(1:end))
xlabel('Geometric mean firing rate')
ylabel({'Synchrony', '(coincidences/spike)'})
title('model ICx')
xlim([0, 60])
