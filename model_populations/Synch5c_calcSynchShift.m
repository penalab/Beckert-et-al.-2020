function Synch5c_calcSynchShift(dir_name, var_filename)

eval(['load ' dir_name 'SynchResults\data\testingvars\' var_filename])

itds = ITD_input;
list = cell2mat(synch_varying(:, 1));

synch_constant = cell(size(list, 1), 7);

tic
for i = 1:size(list, 1)
    synch_constant(i, :) = model_out_synch_calc( ...
        cellfun(@(x) x/1000, spikes2(list(i, 1), :), 'UniformOutput', 0), ...
        cellfun(@(x) x/1000, spikes2(list(i, 2), :), 'UniformOutput', 0), dur/1000);
    synch_constant{i, 1} = list(i, :);
end
toc

win = 5;

syn_con = cellfun(@(x) sum(x{1}(round(length(x{1})/2) - win : round(length(x{1})/2) + win)) / (win*2+1), synch_constant(:, 5:7));

nunits = length(unique(list));
units = nan(nunits);
syn_c = cell(3, 1);

for t = 1:3
    syn_c{t} = nan(nunits);
for i = 1:nunits
    tmp = list == i;
    units(i, :) = [list(tmp(:, 2), 1); nan; list(tmp(:, 1), 2)];
    syn_c{t}(i, :) = [syn_con(tmp(:, 2), t); nan; syn_con(tmp(:, 1), t)];
end

end

eval(['save ' dir_name 'SynchResults\data\out_synch\' var_filename '_syncalc'])

