clear;

% % We are going to shift ICcl spike times to reduce synchrony. Each neuron has a unique constant shift

max_shift = 20; %Time of max shift in ms

% % Simulation parameters

%Stimulus ITD                           %% Can change this value %%
ITD = 100;

if ~exist('repet', 'var')
    repet = 1;
end

r_txt = sprintf('%03d', repet);

% % Directory where SynchResults folder is found

if ~exist('dir_name', 'var')
    dir_name = 'K:\python\model_icx_codeonly\'; %% Can change this value %%
end

eval(['cd ' dir_name '\SynchResults'])

% % Load ICcl spikes from Synch3

eval(['load ' dir_name 'SynchResults\data\iccl\SynchICcl_ITD' num2str(ITD) '_' r_txt ' spikesICcl T BITD CF'])

% This will reduce the number of spikes in ICcl if desired
% spikesICcl = cellfun(@(x) x(1:length(x)/2), spikesICcl, 'UniformOutput', 0);

% This will reduce the stimulus length in ICcl if desired
dur = 300;

[Nf,Nn,Ntrial] = size(spikesICcl);

% % Set up time shifts

Shifts = (rand(Nf,Nn) - 0.5) * max_shift;
% Shifts = (rand(1, Nn) - 0.5) * max_shift;

% % Shift spike times

spikesICcl_shift = spikesICcl;

for n = 1:Nf
    for m = 1:Nn
        for k = 1:Ntrial
            spikesICcl_shift{n,m,k}  = spikesICcl{n,m,k} + Shifts(n,m);
        end
    end
end

spikesICcl_shift = cellfun(@(x) x(x > 0 & x < dur), spikesICcl_shift, 'UniformOutput', 0);

% % Plot rasters of some neurons

figure(1);clf;
for n = 1:9
    subplot(3,3,n);
    bsRaster(squeeze(spikesICcl(n+30,33,:)),[0 T])
    xlim([0 T])
    if n == 2
    title('Example ICcl Rasters')
    end
    xlabel('Time (ms)','fontsize',15)
    ylabel('Trials','fontsize',15)
end
            
% % Plot rasters of same neurons, now shifted

figure(2);clf;
for n = 1:9
    subplot(3,3,n);
    bsRaster(squeeze(spikesICcl_shift(n+30,33,:)),[0 T])
    xlim([0 T])
    if n == 2
    title('Example ICcl Rasters with Shifts')
    end
    xlabel('Time (ms)','fontsize',15)
    ylabel('Trials','fontsize',15)
end

% % Measure synchrony (add code here)


% list = combnk(1:size(spikesICcl_shift, 1), 2);

% list = [1:size(spikesICcl_shift, 1) - 1; 2:size(spikesICcl_shift, 1)]';
i = 1;
list = [ones(size(spikesICcl_shift, 2) - 1, 1) * i, [1:i-1, i+1:size(spikesICcl_shift, 2)]'];

% synch_orig = cell(size(list, 1), 7);
% synch_shift = cell(size(list, 1), 7);
% t = 10;
% for i = 1:size(list, 1)
%     tic
%     disp(['unit 1 = ' num2str(list(i, 1)) '; unit 2 = ' num2str(list(i, 2)) ])
%     synch_orig(i, :) = model_out_synch_calc( ...
%         cellfun(@(x) x/1000, squeeze(spikesICcl(t, list(i, 1), :))', 'UniformOutput', 0), ...
%         cellfun(@(x) x/1000, squeeze(spikesICcl(t, list(i, 2), :))', 'UniformOutput', 0) );
%     synch_shift(i, :) = model_out_synch_calc( ...
%         cellfun(@(x) x/1000, squeeze(spikesICcl_shift(t, list(i, 1), :))', 'UniformOutput', 0), ...
%         cellfun(@(x) x/1000, squeeze(spikesICcl_shift(t, list(i, 2), :))', 'UniformOutput', 0) );
%     synch_orig{i, 1} = list(i, :);
%     synch_shift{i, 1} = list(i, :);
%     toc
% end
  
% load('synch_acrossitd_orig.mat')
synch_acrossitd_orig = cell(size(list, 1), length(CF));

synch_acrossitd_shift = cell(size(list, 1), length(CF));

for t = 1:length(CF)
tic
disp(num2str(CF(t)))
for i = 1:size(list, 1)
    disp(['unit 1 = ' num2str(list(i, 1)) '; unit 2 = ' num2str(list(i, 2)) ])
    st_u1 = cellfun(@(x) x/1000, spikesICcl(t, list(i, 1), :), 'UniformOutput', 0);
    st_u1 = reshape(st_u1, size(st_u1, 2), size(st_u1, 3));
    st_u2 = cellfun(@(x) x/1000, spikesICcl(t, list(i, 2), :), 'UniformOutput', 0);
    st_u2 = reshape(st_u2, size(st_u2, 2), size(st_u2, 3));
    synch_acrossitd_orig{i, t} = model_out_synch_calc(st_u1, st_u2, dur/1000);
    synch_acrossitd_orig{i, t}{1} = list(i, :);
    st_u1 = cellfun(@(x) x/1000, spikesICcl_shift(t, list(i, 1), :), 'UniformOutput', 0);
    st_u1 = reshape(st_u1, size(st_u1, 2), size(st_u1, 3));
    st_u2 = cellfun(@(x) x/1000, spikesICcl_shift(t, list(i, 2), :), 'UniformOutput', 0);
    st_u2 = reshape(st_u2, size(st_u2, 2), size(st_u2, 3));
    synch_acrossitd_shift{i, t} = model_out_synch_calc(st_u1, st_u2, dur/1000);
    synch_acrossitd_shift{i, t}{1} = list(i, :);
end
toc
end

% %

win = 1;

gm_acrossitd_orig = cellfun(@(x) sqrt(x{2} .* x{3}), synch_acrossitd_orig);
ccg_acrossitd_orig = cellfun(@(x) sum(x{5}{1}(round(length(x{5}{1})/2) - win : round(length(x{5}{1})/2) + win)) / (win*2+1), synch_acrossitd_orig);
sh_acrossitd_orig = cellfun(@(x) sum(x{6}{1}(round(length(x{6}{1})/2) - win : round(length(x{6}{1})/2) + win)) / (win*2+1), synch_acrossitd_orig);
cor_acrossitd_orig = cellfun(@(x) sum(x{7}{1}(round(length(x{7}{1})/2) - win : round(length(x{7}{1})/2) + win)) / (win*2+1), synch_acrossitd_orig);
list_acrossitd_orig = cell2mat(cellfun(@(x) x{1}, synch_acrossitd_orig(:, 1), 'UniformOutput', 0));

gm_acrossitd_shift = cellfun(@(x) sqrt(x{2} .* x{3}), synch_acrossitd_shift);
ccg_acrossitd_shift = cellfun(@(x) sum(x{5}{1}(round(length(x{5}{1})/2) - win : round(length(x{5}{1})/2) + win)) / (win*2+1), synch_acrossitd_shift);
sh_acrossitd_shift = cellfun(@(x) sum(x{6}{1}(round(length(x{6}{1})/2) - win : round(length(x{6}{1})/2) + win)) / (win*2+1), synch_acrossitd_shift);
cor_acrossitd_shift = cellfun(@(x) sum(x{7}{1}(round(length(x{7}{1})/2) - win : round(length(x{7}{1})/2) + win)) / (win*2+1), synch_acrossitd_shift);
list_acrossitd_shift = cell2mat(cellfun(@(x) x{1}, synch_acrossitd_shift(:, 1), 'UniformOutput', 0));

% %

sim = nan(10000, 4);
dif = nan(10000, 4);

for f = 1:size(ccg_acrossitd_shift, 2)
    for i = 1:size(ccg_acrossitd_shift, 1)
        if abs(list_acrossitd_shift(i, 1) - list_acrossitd_shift(i, 2)) < 3
            L = find(isnan(sim(:, 1)), 1, 'first');
            sim(L, 1) = gm_acrossitd_shift(i, f);
            sim(L, 2) = ccg_acrossitd_shift(i, f);
            sim(L, 3) = sh_acrossitd_shift(i, f);
            sim(L, 4) = cor_acrossitd_shift(i, f);
        elseif abs(list_acrossitd_shift(i, 1) - list_acrossitd_shift(i, 2)) > 50
            L = find(isnan(dif(:, 1)), 1, 'first');
            dif(L, 1) = gm_acrossitd_shift(i, f);
            dif(L, 2) = ccg_acrossitd_shift(i, f);
            dif(L, 3) = sh_acrossitd_shift(i, f);
            dif(L, 4) = cor_acrossitd_shift(i, f);
        end
        
    end
end

L = find(isnan(sim(:, 1)), 1, 'first');
sim(L:end, :) = [];
L = find(isnan(dif(:, 1)), 1, 'first');
dif(L:end, :) = [];

%%

load('model_synch_iccl_shifttest_10_shortdur.mat')

figure
set(gcf, 'position', [0 0 1000 300])

subplot(1, 10, 1:7)
boxplot2_points(1:2, {sim(:, 2), dif(:, 2)}, 0.75)
boxplot2_points(4:5, {sim(:, 3), dif(:, 3)}, 0.75)
boxplot2_points(7:8, {sim(:, 4), dif(:, 4)}, 0.75)
xticks([1, 2, 4, 5, 7, 8])
xlim([0, 9])
ylabel('Synchrony (coincidences/spikes)')

% %

a = gm_acrossitd_shift(:);
colors = {'k', 'r', 'b'};
subplot(1, 10, 8:10)
hold on

for t = 1:3
    switch t
        case 1
            b = ccg_acrossitd_shift(:);
        case 2
            b = sh_acrossitd_shift(:);
        case 3
            b = cor_acrossitd_shift(:);
    end

scatter(a(1:end), b(1:end), 10, 'fill', colors{t})
xlabel('Geometric mean firing rate')
ylabel('Synchrony (coincidences/spike)')

% [cc, pp] = corr(a,b,'rows','pairwise');
% title(['r = ' num2str(cc) '; p = ' num2str(pp)])

end


