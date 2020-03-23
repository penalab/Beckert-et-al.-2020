load('f:\desktop\Code\Reproducibility\python\reproducibility\data\ICxSpikes_ITD0.mat')

subplot(1,2,1)

load('K:\python\reproducibility\data\repro_iccl_ITD100_002_reducedSNR.mat')
load('K:\python\reproducibility\data\repro_iccl_ITD100_002_reducedSNR_strength.mat')

pa = unique(params);

repro = nan(length(pa), 1);

for p = 1:length(pa)
    repro(p, 1) = nanmean(strength(params == pa(p)));
end

imagesc(reshape(repro, 5, 3))
xticks(1:3)
xticklabels(unique(ge))
xlabel('ge')
yticks(1:5)
yticklabels(unique(gi))
ylabel('gi')
h = colorbar;
set(get(h, 'label'), 'string', 'reproducibility')
title('ITD0')
axis square

subplot(1,2,2)

load('f:\desktop\Code\Reproducibility\python\reproducibility\data\repro_ITD150.mat')
load('f:\desktop\Code\Reproducibility\python\reproducibility\data\repro_ITD150strength.mat')

pa = unique(params);

repro = nan(length(pa), 1);

for p = 1:length(pa)
    repro(p, 1) = nanmean(strength(params == pa(p)));
end

imagesc(reshape(repro, 5, 3))
xticks(1:3)
xticklabels(unique(ge))
xlabel('ge')
yticks(1:5)
yticklabels(unique(gi))
ylabel('gi')
h = colorbar;
set(get(h, 'label'), 'string', 'reproducibility')
title('ITD150')
axis square

%%

load('f:\desktop\Code\python\reproducibility\data\ICxSpikes_ITD0.mat')
load('f:\desktop\Code\python\reproducibility\data\repro_ITD0.mat')
load('f:\desktop\Code\python\reproducibility\data\repro_ITD0strength.mat')

for c = 1:size(curve, 1)
    curve(c, :) = (curve(c, :) - fr(c) ^ 2) * cw / fr(c);
end

center = round(size(curve, 2) / 2);
win = (size(curve, 2) - center) * cw;
x = -win:cw:win;

figure(1)
plot(x, curve)

cur = cell(length(unique(params)), 1);
for p = unique(params)
    cur{p+1} = curve(params == p, :);
end

clear p win center

c = cellfun(@(x) nanmean(x, 1), cur, 'UniformOutput', 0);
c = reshape(c, 5, 3)';
ge = reshape(ge, 5, 3)';
gi = reshape(gi, 5, 3)';

figure(2)
for i = 1:length(cur)
    subplot(5, 3, i)
    plot(x, c{i})
    xlabel(num2str(ge(i)))
    ylabel(num2str(gi(i)))
    xlim([-0.02, 0.02])
end


