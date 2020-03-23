itds = [-180:20:-80, 80:20:180];

dir_name = 'K:\python\model_icx_codeonly\';

load([dir_name 'SynchResults\data\model_bounds'], 'bounds', 'itd');

b_out = nan(length(itds), 4);

mani = cell(length(itds), 100);
orig = cell(length(itds), 100);

for t = 1:length(itds)
    
    ix = find(itd == itds(t), 1, 'first');
    b_out(t, :) = bounds(ix, :);
    
    load([dir_name 'SynchResults\data\out_mani\ManipulateICx_ITD' num2str(itds(t))])
    
    for r = 1:size(sc, 1)
        mani{t, r} = sc2(r, :);
        orig{t, r} = sc(r, :);
    end
    
end

output_mani = mani;
output_orig = orig;
itd = itds;
bounds = b_out;

save model_out_full_mani output_orig output_mani itd bounds

%% Run model_out_synch_strength

itds = [-180:20:-80, 80:20:180];

dir_name = 'K:\python\model_icx_codeonly\';

rep = 1;

for ITD = itds
    model_out_synch_strength(ITD, rep, dir_name);
end
