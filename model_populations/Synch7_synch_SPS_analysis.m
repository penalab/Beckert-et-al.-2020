cd('K:\python\model_icx_codeonly\SynchResults\data\out_synch')

files = ls;
files = files(3:end,:);

syn_curve_orig_ccg = nan(size(files, 1), 81);
syn_curve_mani_ccg = nan(size(files, 1), 81);
syn_curve_orig_sh = nan(size(files, 1), 81);
syn_curve_mani_sh = nan(size(files, 1), 81);

sidepeaks = nan(size(files, 1), 2);

for f = 1:size(files, 1)
    
    load(files(f, :), 'SC', 'SC2', 'syn_v', 'syn_c', 'ind2')
    
    M = mean(SC,1);
    M2 = mean(SC2,1);
    sidepeaks(f, 1) = max(M(ind2))/max(M);
    sidepeaks(f, 2) = max(M2(ind2))/max(M2);
    
    syn_curve_orig_ccg(f, 1) = nanmean(syn_v{1}(1, 1:2));
    syn_curve_mani_ccg(f, 1) = nanmean(syn_c{1}(1, 1:2));
    syn_curve_orig_ccg(f, end) = nanmean(syn_v{1}(end, end-1:end));
    syn_curve_mani_ccg(f, end) = nanmean(syn_c{1}(end, end-1:end));
    
    for p = 2:80
        syn_curve_orig_ccg(f, p) = nanmean(syn_v{1}(p, p-1:p+1));
        syn_curve_mani_ccg(f, p) = nanmean(syn_c{1}(p, p-1:p+1));
    end
    
    syn_curve_orig_sh(f, 1) = nanmean(syn_v{2}(1, 1:2));
    syn_curve_mani_sh(f, 1) = nanmean(syn_c{2}(1, 1:2));
    syn_curve_orig_sh(f, end) = nanmean(syn_v{2}(end, end-1:end));
    syn_curve_mani_sh(f, end) = nanmean(syn_c{2}(end, end-1:end));
    
    for p = 2:80
        syn_curve_orig_sh(f, p) = nanmean(syn_v{2}(p, p-1:p+1));
        syn_curve_mani_sh(f, p) = nanmean(syn_c{2}(p, p-1:p+1));
    end  
    
end

clear M M2 f p SC SC2 syn_v syn_c ind2


%%

itd = linspace(-200, 200, 81);

for f = 1:size(sidepeaks, 1)
    
    figure(f)
    
    subplot(1, 2, 1)
    plot(itd, syn_curve_orig_ccg(f, :), 'b')
    hold on
    plot(itd, syn_curve_mani_ccg(f, :), 'r')
    
    subplot(1, 2, 2)
    plot(itd, syn_curve_orig_sh(f, :), 'b')
    hold on
    plot(itd, syn_curve_mani_sh(f , :), 'r')
    
end





