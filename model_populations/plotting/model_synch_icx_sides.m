itds = -200:20:200;

rep = 1;
r_txt = sprintf('%03d', rep);
% dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\';
dir_name = 'K:\python\model_icx_codeonly\';
load([dir_name 'SynchResults\model_out_full_mani.mat'])

peaks_s = bounds(:, 3) + round((bounds(:, 4) - bounds(:, 3)) / 2);
peaks_m = bounds(:, 1) + round((bounds(:, 2) - bounds(:, 1)) / 2);

main = cell(length(itd), 4);
side_c = cell(length(itd), 4);
side_v = cell(length(itd), 4);

C_main = cell(length(itd), 4);
C_side = cell(length(itd), 4);
V_main = cell(length(itd), 4);
V_side = cell(length(itd), 4);

for ITD = 1:length(itd)
    
    load([dir_name 'SynchResults\data\out_mani\SynchOut_strength_ITD' num2str(itd(ITD)) '_' r_txt])
    
    g = nanmean([[nan; diag(gm, -1)],[diag(gm, 1); nan]], 2);
    C_main{ITD, 1} = g(bounds(ITD, 1): bounds(ITD, 2));
    C_side{ITD, 1} = g(bounds(ITD, 3): bounds(ITD, 4));
    V_main{ITD, 1} = g(bounds(ITD, 1): bounds(ITD, 2));
    V_side{ITD, 1} = g(bounds(ITD, 3): bounds(ITD, 4));
    
    main{ITD, 1} = g(peaks_m(ITD)-1:peaks_m(ITD)+1);
    if peaks_s(ITD) ~= 1
        side_c{ITD, 1} = g(peaks_s(ITD)-1:peaks_s(ITD)+1);
        side_v{ITD, 1} = g(peaks_s(ITD)-1:peaks_s(ITD)+1);
    else
        side_c{ITD, 1} = nan;
        side_v{ITD, 1} = nan;
    end
    
    for t = 1:3
        
        c = nanmean([[nan; diag(syn_c{t}, -1)],[diag(syn_c{t}, 1); nan]], 2);
        v = nanmean([[nan; diag(syn_v{t}, -1)],[diag(syn_v{t}, 1); nan]], 2);
        
        C_main{ITD, t+1} = c(bounds(ITD, 1): bounds(ITD, 2));
        C_side{ITD, t+1} = c(bounds(ITD, 3): bounds(ITD, 4));
        V_main{ITD, t+1} = v(bounds(ITD, 1): bounds(ITD, 2));
        V_side{ITD, t+1} = v(bounds(ITD, 3): bounds(ITD, 4));
        
        main{ITD, t+1} = c(peaks_m(ITD)-1:peaks_m(ITD)+1);
        if peaks_s(ITD) ~= 1
            side_c{ITD, t+1} = c(peaks_s(ITD)-1:peaks_s(ITD)+1);
            side_v{ITD, t+1} = v(peaks_s(ITD)-1:peaks_s(ITD)+1);
        else
            side_c{ITD, t+1} = nan;
            side_v{ITD, t+1} = nan;
        end
        
    end
    
end

C_main = cell2mat(C_main);
C_side = cell2mat(C_side);
V_main = cell2mat(V_main);
V_side = cell2mat(V_side);
main = cell2mat(main);
side_c = cell2mat(side_c);
side_v = cell2mat(side_v);

%%

names = {'standard', 'shifted', 'corrected'};
figure(3)
set(gcf, 'Position', [0 0 1200 800])

for t = 2:4
    subplot(2, 3, t-1)
    scatter(C_main(:, 1), C_main(:, t), 10, 'k', 'fill')
    hold on
    scatter(C_side(:, 1), C_side(:, t), 10, 'r', 'fill')
    scatter(V_side(:, 1), V_side(:, t), 10, 'b', 'fill')
    
    xlabel('Geometric mean firing rate')
    ylabel('Synchrony (coincidences/spike)')
    
    m = corr(C_main(:, 1), C_main(:, t), 'rows', 'pairwise');
    c = corr(C_side(:, 1), C_side(:, t), 'rows', 'pairwise');
    v = corr(V_side(:, 1), V_side(:, t), 'rows', 'pairwise');
    
    title({names{t-1}, ['main peak = ' num2str(m)], ...
        ['orig side peak = ' num2str(v)], ...
        ['manipulated side peak = ' num2str(c)]});

end


for t = 2:4
    subplot(2, 3, t+2)
    boxplot2_points(1:3, {main(:, t), side_v(:, t), side_c(:, t)}, 0.75);
    xlim([0 4])
    xticks(1:3)
    xticklabels({'Main peak', 'Original side peak', 'Manipulated side peak'})
    xtickangle(45)
    ylabel('Synchrony (coincidences/spike)')
    title(names{t-1})
end

