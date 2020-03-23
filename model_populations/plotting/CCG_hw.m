load('K:\WorkingFolder\synchony\ot_synch_ff_del_fr_200_cor.mat')

gm = ot_synch.ff{3} .* ot_synch.ff{4};
synch = cellfun(@(x) max(x(100)), ot_synch.ff{11});

%%

xx = -100:1:100;
xx = xx(2:end-1);

hw = nan(size(gm, 1), 1);
curves = nan(size(gm, 1), length(xx));

for i = 1:size(gm, 1)
    [~, ix] = max(synch(i, :));
    curves(i, :) = ot_synch.ff{11}{i, ix}/max(ot_synch.ff{11}{i, ix});
    if max(curves(i, :)) > 0
        half = max(curves(i, :))/2;
        [~, ix] = min(abs(curves(i, :) - half));
        hw(i) = abs(xx(ix));
    end
end

clear i ix curve half xx

imagesc(curves)

%%

hw(hw == 99) = [];
hw(isnan(hw)) = [];

median(hw)