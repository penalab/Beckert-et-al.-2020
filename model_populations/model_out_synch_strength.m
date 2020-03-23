function model_out_synch_strength(ITD, rep, dir_name)
%%
if ~exist('ITD', 'var')
    ITD = 100;
end
if ~exist('rep', 'var') 
    rep = 1;
end
r_txt = sprintf('%03d', rep);

if ~exist('dir_name', 'var')
    dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\';
end

load([dir_name 'SynchResults\data\SynchParms'], 'bestITD')
load([dir_name 'SynchResults\data\out_old\SynchOut_ITD' num2str(ITD) '_' r_txt])
% load([dir_name 'SynchResults\data\tmp\SynchOut_ITD' num2str(ITD) '_' r_txt])
load([dir_name 'SynchResults\data\out_mani\ManipulateICx_ITD' num2str(ITD)])
%%

list = combnk(1:size(spikesICx, 1), 2);

synch_varying = cell(size(list, 1), 7);
synch_constant = cell(size(list, 1), 7);

tic
for i = 1:size(list, 1)
    disp(['unit 1 = ' num2str(list(i, 1)) '; unit 2 = ' num2str(list(i, 2)) ])
    synch_varying(i, :) = model_out_synch_calc( ...
        cellfun(@(x) x/1000, spikesICx(list(i, 1), :), 'UniformOutput', 0), ...
        cellfun(@(x) x/1000, spikesICx(list(i, 2), :), 'UniformOutput', 0) );
    synch_constant(i, :) = model_out_synch_calc( ...
        cellfun(@(x) x/1000, spikes2(list(i, 1), :), 'UniformOutput', 0), ...
        cellfun(@(x) x/1000, spikes2(list(i, 2), :), 'UniformOutput', 0) );
    synch_varying{i, 1} = list(i, :);
    synch_constant{i, 1} = list(i, :);
end
toc

% %

% win = 10;
win = 1;

list = cell2mat(synch_varying(:, 1));

gm_con = sqrt(cell2mat(synch_constant(:, 2)) .* cell2mat(synch_constant(:, 3)));
gm_var = sqrt(cell2mat(synch_varying(:, 2)) .* cell2mat(synch_varying(:, 3)));

syn_con = cellfun(@(x) sum(x{1}(round(length(x{1})/2) - win : round(length(x{1})/2) + win)) / (win*2+1), synch_constant(:, 5:7));
syn_var = cellfun(@(x) sum(x{1}(round(length(x{1})/2) - win : round(length(x{1})/2) + win)) / (win*2+1), synch_varying(:, 5:7));

% %

nunits = length(unique(list));
units = nan(nunits);
syn_c = cell(3, 1);
syn_v = cell(3, 1);

for t = 1:3
    syn_c{t} = nan(nunits);
    syn_v{t} = nan(nunits);
    gm = nan(nunits);
for i = 1:nunits
    tmp = list == i;
    units(i, :) = [list(tmp(:, 2), 1); nan; list(tmp(:, 1), 2)];
    syn_c{t}(i, :) = [syn_con(tmp(:, 2), t); nan; syn_con(tmp(:, 1), t)];
    syn_v{t}(i, :) = [syn_var(tmp(:, 2), t); nan; syn_var(tmp(:, 1), t)];
    gm(i, :) = [gm_var(tmp(:, 2), 1); nan; gm_var(tmp(:, 1), 1)];
end

end

clear i t

%%

save([dir_name 'SynchResults\data\out_mani\SynchOut_strength_ITD' num2str(ITD) '_' r_txt], ...
    'spikesICx', 'spikes2', 'gm', 'syn_c', 'syn_v', 'gm_con', 'gm_var', 'syn_con', ...
    'syn_var', 'bestITD', 'synch_constant', 'synch_varying', 'win', 'nunits')
% save([dir_name 'SynchResults\data\tmp\SynchOut_strength_ITD' num2str(ITD) '_' r_txt], ...
%     'spikesICx', 'spikes2', 'output_orig', 'output_mani', 'gm', 'syn_c', 'syn_v', 'gm_con', 'gm_var', 'syn_con', ...
%     'syn_var', 'bestITD', 'synch_constant', 'synch_varying', 'win', 'nunits')
end