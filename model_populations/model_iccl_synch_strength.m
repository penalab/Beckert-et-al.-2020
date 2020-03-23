function model_iccl_synch_strength(ITD, rep, dir_name)

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

load([dir_name 'SynchResults\data\iccl\SynchICcl_ITD' num2str(ITD) '_' r_txt '_SNR20'], 'BITD', 'CF', 'spikesICcl')

%%

% list = combnk(1:size(spikesICcl, 1), 2);
% synch_acrossfreq = cell(size(list, 1), length(BITD));
% 
% for t = 1:length(BITD)
% tic
% disp(num2str(BITD(t)))
% for i = 1:size(list, 1)
%     disp(['unit 1 = ' num2str(list(i, 1)) '; unit 2 = ' num2str(list(i, 2)) ])
%     st_u1 = cellfun(@(x) x/1000, spikesICcl(list(i, 1), t, :), 'UniformOutput', 0);
%     st_u1 = reshape(st_u1, size(st_u1, 2), size(st_u1, 3));
%     st_u2 = cellfun(@(x) x/1000, spikesICcl(list(i, 2), t, :), 'UniformOutput', 0);
%     st_u2 = reshape(st_u2, size(st_u2, 2), size(st_u2, 3));
%     synch_acrossfreq{i, t} = model_out_synch_calc(st_u1, st_u2);
%     synch_acrossfreq{i, t}{1} = list(i, :);
% end
% toc
% end


list = combnk(1:size(spikesICcl, 2), 2);
synch_acrossitd = cell(size(list, 1), length(CF));

for t = 1:length(CF)
tic
disp(num2str(CF(t)))
for i = 1:size(list, 1)
    disp(['unit 1 = ' num2str(list(i, 1)) '; unit 2 = ' num2str(list(i, 2)) ])
    st_u1 = cellfun(@(x) x/1000, spikesICcl(t, list(i, 1), :), 'UniformOutput', 0);
    st_u1 = reshape(st_u1, size(st_u1, 2), size(st_u1, 3));
    st_u2 = cellfun(@(x) x/1000, spikesICcl(t, list(i, 2), :), 'UniformOutput', 0);
    st_u2 = reshape(st_u2, size(st_u2, 2), size(st_u2, 3));
    synch_acrossitd{i, t} = model_out_synch_calc(st_u1, st_u2);
    synch_acrossitd{i, t}{1} = list(i, :);
end
toc
end

clear t i list st_u1 st_u2
%%

win = 5;

% gm_acrossfreq = cellfun(@(x) sqrt(x{2} .* x{3}), synch_acrossfreq);
gm_acrossitd = cellfun(@(x) sqrt(x{2} .* x{3}), synch_acrossitd);

% ccg_acrossfreq = cellfun(@(x) sum(x{5}{1}(round(length(x{5}{1})/2) - win : round(length(x{5}{1})/2) + win)) / (win*2+1), synch_acrossfreq);
% sh_acrossfreq = cellfun(@(x) sum(x{6}{1}(round(length(x{6}{1})/2) - win : round(length(x{6}{1})/2) + win)) / (win*2+1), synch_acrossfreq);
% cor_acrossfreq = cellfun(@(x) sum(x{7}{1}(round(length(x{7}{1})/2) - win : round(length(x{7}{1})/2) + win)) / (win*2+1), synch_acrossfreq);

ccg_acrossitd = cellfun(@(x) sum(x{5}{1}(round(length(x{5}{1})/2) - win : round(length(x{5}{1})/2) + win)) / (win*2+1), synch_acrossitd);
sh_acrossitd = cellfun(@(x) sum(x{6}{1}(round(length(x{6}{1})/2) - win : round(length(x{6}{1})/2) + win)) / (win*2+1), synch_acrossitd);
cor_acrossitd = cellfun(@(x) sum(x{7}{1}(round(length(x{7}{1})/2) - win : round(length(x{7}{1})/2) + win)) / (win*2+1), synch_acrossitd);

% list_acrossfreq = cell2mat(cellfun(@(x) x{1}, synch_acrossfreq(:, 1), 'UniformOutput', 0));
list_acrossitd = cell2mat(cellfun(@(x) x{1}, synch_acrossitd(:, 1), 'UniformOutput', 0));

%%

% save([dir_name 'SynchResults\data\iccl_strength\SynchICcl_strength_ITD' num2str(ITD) '_' r_txt], ...
%     'synch_acrossfreq', 'ccg_acrossfreq', 'sh_acrossfreq', 'cor_acrossfreq', 'list_acrossfreq', 'gm_acrossfreq', ...
%     'synch_acrossitd', 'syn_acrossitd', 'sh_acrossitd', 'cor_acrossitd', 'list_acrossitd', 'gm_acrossitd');
save([dir_name 'SynchResults\data\iccl_strength\SynchICcl_strength_ITD' num2str(ITD) '_' r_txt '_SNR20_itd'], ...
    'synch_acrossitd', 'ccg_acrossitd', 'sh_acrossitd', 'cor_acrossitd', 'list_acrossitd', 'gm_acrossitd');
% save([dir_name 'SynchResults\data\iccl_strength\SynchICcl_strength_ITD' num2str(ITD) '_' r_txt '_freq'], ...
%     'synch_acrossfreq', 'ccg_acrossfreq', 'sh_acrossfreq', 'cor_acrossfreq', 'list_acrossfreq', 'gm_acrossfreq');
end