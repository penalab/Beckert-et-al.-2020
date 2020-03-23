cd('K:\WorkingFolder')
pathstartup_pc
cd('K:\python\model_icx_codeonly\SynchResults')
load('K:\python\model_icx_codeonly\SynchResults\data\SynchICcl_Filt_ITD0.mat')

%%

[Nfreq, Nneu, Nfilt, Nrep] = size(spikesICcl);

list = combnk(1:5:Nneu, 2);
synch_orig = cell(size(list, 1), Nfreq);
synch_filt = cell(size(list, 1), Nfreq);

for t = 1:length(CF)
tic
disp(num2str(CF(t)))
for i = 1:size(list, 1)
    disp(['unit 1 = ' num2str(list(i, 1)) '; unit 2 = ' num2str(list(i, 2)) ])
    st_u1 = cellfun(@(x) x/1000, spikesICcl(t, list(i, 1), 1, :), 'UniformOutput', 0);
    st_u1 = reshape(st_u1, size(st_u1, 2), size(st_u1, 4));
    st_u2 = cellfun(@(x) x/1000, spikesICcl(t, list(i, 2), 1, :), 'UniformOutput', 0);
    st_u2 = reshape(st_u2, size(st_u2, 2), size(st_u2, 4));
    synch_orig{i, t} = model_out_synch_calc(st_u1, st_u2);
    synch_orig{i, t}{1} = list(i, :);
    st_u1 = cellfun(@(x) x/1000, spikesICcl(t, list(i, 1), 1, :), 'UniformOutput', 0);
    st_u1 = reshape(st_u1, size(st_u1, 2), size(st_u1, 4));
    st_u2 = cellfun(@(x) x/1000, spikesICcl(t, list(i, 2), 3, :), 'UniformOutput', 0);
    st_u2 = reshape(st_u2, size(st_u2, 2), size(st_u2, 4));
    synch_filt{i, t} = model_out_synch_calc(st_u1, st_u2);
    synch_filt{i, t}{1} = list(i, :);
end
toc
end

clear t i list st_u1 st_u2
%%

win = 10;

gm = cellfun(@(x) sqrt(x{2} .* x{3}), synch_orig);

ccg_orig = cellfun(@(x) sum(x{5}{1}(round(length(x{5}{1})/2) - win : round(length(x{5}{1})/2) + win)) / (win*2+1), synch_orig);
sh_orig = cellfun(@(x) sum(x{6}{1}(round(length(x{6}{1})/2) - win : round(length(x{6}{1})/2) + win)) / (win*2+1), synch_orig);
cor_orig = cellfun(@(x) sum(x{7}{1}(round(length(x{7}{1})/2) - win : round(length(x{7}{1})/2) + win)) / (win*2+1), synch_orig);

ccg_filt = cellfun(@(x) sum(x{5}{1}(round(length(x{5}{1})/2) - win : round(length(x{5}{1})/2) + win)) / (win*2+1), synch_filt);
sh_filt = cellfun(@(x) sum(x{6}{1}(round(length(x{6}{1})/2) - win : round(length(x{6}{1})/2) + win)) / (win*2+1), synch_filt);
cor_filt = cellfun(@(x) sum(x{7}{1}(round(length(x{7}{1})/2) - win : round(length(x{7}{1})/2) + win)) / (win*2+1), synch_filt);


