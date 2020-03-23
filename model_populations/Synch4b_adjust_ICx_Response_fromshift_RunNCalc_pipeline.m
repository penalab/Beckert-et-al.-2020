function Synch4b_adjust_ICx_Response_fromshift_RunNCalc_pipeline(repet, ITD_input, dir_name, SNR_ICcl, SNR_ICx, ICx_weight, ICx_side, ICx_width, freq_width, BW, gainL, gainU, W_weight)

tic

if ~exist('repet', 'var')
    repet = 0;
end

r_txt = sprintf('%03d', repet);

%% Simulation parameters

%Stimulus ITD                           %% Can change this value %%
if ~exist('ITD_input', 'var')
    ITD_input = 0;
end

eval(['load ' dir_name 'SynchResults\data\icx\SynchICx_ITD' num2str(ITD_input) '_' r_txt ...
    '_SNRICcl' num2str(SNR_ICcl * 10) ...
    '_SNRICx' num2str(SNR_ICx * 100000) ...
    '_ICxweight' num2str(ICx_weight * 10) ...
    '_ICxside' num2str(ICx_side * 10) ...
    '_ICxwidth' num2str(ICx_width * 10) ...
    '_freqwidth' num2str(freq_width * 10) ...
    '_bw' num2str(BW * 100000) ...
    '_gainL' num2str(gainL * 100) ...
    '_gainU' num2str(gainU)])

% shuffle spike trains

for i = 1:size(spikesICx, 1)
    spikesICx(i, :) = spikesICx(i, randperm(size(spikesICx, 2)));
end

max_shift = 0;
dur = 100;  %%changed (originally 300)

[Nn, Ntrial] = size(spikesICx);
Shifts = (rand(Nn,1) - 0.5) * max_shift;

% % Shift spike times

for m = 1:Nn
    for k = 1:Ntrial
        spikesICx{m,k}  = spikesICx{m,k} + Shifts(m,1);
    end
end
spikesICx = cellfun(@(x) x(x > 150 & x < dur+150), spikesICx, 'UniformOutput', 0);

itds = ITD_input;

list = combnk(1:size(spikesICx, 1), 2);

synch_varying = cell(size(list, 1), 7);
% synch_constant = cell(size(list, 1), 7);

tic
for i = 1:size(list, 1)
    synch_varying(i, :) = model_out_synch_calc( ...
        cellfun(@(x) x/1000, spikesICx(list(i, 1), :), 'UniformOutput', 0), ...
        cellfun(@(x) x/1000, spikesICx(list(i, 2), :), 'UniformOutput', 0), dur/1000);
%     synch_constant(i, :) = model_out_synch_calc( ...
%         cellfun(@(x) x/1000, spikes2(list(i, 1), :), 'UniformOutput', 0), ...
%         cellfun(@(x) x/1000, spikes2(list(i, 2), :), 'UniformOutput', 0) );
    synch_varying{i, 1} = list(i, :);
%     synch_constant{i, 1} = list(i, :);
end
toc

% %

win = 5;

list = cell2mat(synch_varying(:, 1));

% gm_con = sqrt(cell2mat(synch_constant(:, 2)) .* cell2mat(synch_constant(:, 3)));
gm_var = sqrt(cell2mat(synch_varying(:, 2)) .* cell2mat(synch_varying(:, 3)));

% syn_con = cellfun(@(x) sum(x{1}(round(length(x{1})/2) - win : round(length(x{1})/2) + win)) / (win*2+1), synch_constant(:, 5:7));
syn_var = cellfun(@(x) sum(x{1}(round(length(x{1})/2) - win : round(length(x{1})/2) + win)) / (win*2+1), synch_varying(:, 5:7));

% %

nunits = length(unique(list));
units = nan(nunits);
% syn_c = cell(3, 1);
syn_v = cell(3, 1);

for t = 1:3
%     syn_c{t} = nan(nunits);
    syn_v{t} = nan(nunits);
    gm = nan(nunits);
for i = 1:nunits
    tmp = list == i;
    units(i, :) = [list(tmp(:, 2), 1); nan; list(tmp(:, 1), 2)];
%     syn_c{t}(i, :) = [syn_con(tmp(:, 2), t); nan; syn_con(tmp(:, 1), t)];
    syn_v{t}(i, :) = [syn_var(tmp(:, 2), t); nan; syn_var(tmp(:, 1), t)];
    gm(i, :) = [gm_var(tmp(:, 2), 1); nan; gm_var(tmp(:, 1), 1)];
end

end

clear i t

% %
sim = nan(10000, 4);
dif = nan(10000, 4);

GM = cell(length(itds), 1);
SYN = cell(length(itds), 1);

for ITD = 1:length(itds)

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

% %

figure
set(gcf, 'position', [0 0 1000 300])

subplot(1, 10, 1:7)
boxplot2_points(1:2, {sim(:, 2), dif(:, 2)}, 0.75)
boxplot2_points(4:5, {sim(:, 3), dif(:, 3)}, 0.75)
boxplot2_points(7:8, {sim(:, 4), dif(:, 4)}, 0.75)
xticks([1, 2, 4, 5, 7, 8])
xlim([0, 9])
ylabel({'Synchrony', '(coincidences/spikes)'})


a = cell2mat(GM);
colors = {'k', 'r', 'b'};
syn = cell2mat(SYN);

subplot(1, 10, 8:10)
hold on

for t = 1:3

    b = syn(:, t);

    scatter(a(1:end), b(1:end), 10, 'fill', colors{t})
    xlabel('Geometric mean firing rate')
    ylabel({'Synchrony', '(coincidences/spike)'})

end

eval(['save ' dir_name 'SynchResults\data\icx\SynchICx_ITD' num2str(ITD_input) '_' r_txt ...
    '_SNRICcl' num2str(SNR_ICcl * 10) ...
    '_SNRICx' num2str(SNR_ICx * 100000) ...
    '_ICxweight' num2str(ICx_weight * 10) ...
    '_ICxside' num2str(ICx_side * 10) ...
    '_ICxwidth' num2str(ICx_width * 10) ...
    '_freqwidth' num2str(freq_width * 10) ...
    '_bw' num2str(BW * 100000) ...
    '_gainL' num2str(gainL * 100) ...
    '_gainU' num2str(gainU) ...
    '_shiftshuff'])

toc
