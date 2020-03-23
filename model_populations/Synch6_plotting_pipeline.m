cd('K:\python\model_icx_codeonly\SynchResults\data\testingvars_firstpass')

files = ls;
files = files(3:end,:);
params = cell(size(files, 1), 1);

for f = 1:size(files, 1)
    files(f, regexp(files(f,:),'[modeloutITDSNRICclSNRICxICxweightICxsideICxwidthfreqwidthshiftshuffbw.mat ]')) = '_';
    files(f, regexp(files(f,:),'[AM]')) = '1';
    params{f} = cellfun(@(x) str2double(x)/10, strsplit(files(f, :), '_'));
    if length(params{f}) == 11
        params{f}(end) = [];
        params{f}(end) = 1;
    end
    params{f}(1) = [];
    params{f}(1) = params{f}(1) * 10;
    params{f}(2) = params{f}(2) * 10;
end

params = cell2mat(params);

clear files f

%%

dir_name = 'K:\python\model_icx_codeonly\';  %% Can change this value %%

% Full paramaters tested
itds = unique(params(:, 1))';
repet = unique(params(:, 2))';

snr_iccl = unique(params(:, 3))';
snr_icx = unique(params(:, 4))';
icx_weight = unique(params(:, 5))';
icx_side = unique(params(:, 6))';
icx_width = unique(params(:, 7))';
freq_width = unique(params(:, 8))';

count_total = 0;
count_found = 0;

params_used = nan(500, 9);

for SNR_ICcl = snr_iccl
for SNR_ICx = snr_icx
for ICx_weight = icx_weight
    
for ICx_side = icx_side
for ICx_width = icx_width
for FREQ_width = freq_width
    
    
for rep = repet
for ITD = itds
    
    count_total = count_total + 1;
    
    r_txt = sprintf('%03d', rep);

    if isfile([dir_name 'SynchResults\data\testingvars_firstpass\modelout_ITD' num2str(ITD) '_' r_txt ...
    '_SNRICcl' num2str(SNR_ICcl * 10) ...
    '_SNRICx' num2str(SNR_ICx * 10) ...
    '_ICxweight' num2str(ICx_weight * 10) ...
    '_ICxside' num2str(ICx_side * 10) ...
    '_ICxwidth' num2str(ICx_width * 10) ...
    '_freqwidth' num2str(FREQ_width * 10) ...
    '_shiftshuff.mat'])


    count_found = count_found + 1;

    L = find(isnan(params_used(:, 1)), 1, 'first');
    params_used(L, 1) = ITD;
    params_used(L, 2) = rep;
    params_used(L, 3) = SNR_ICcl;
    params_used(L, 4) = SNR_ICx;
    params_used(L, 5) = ICx_weight;
    params_used(L, 6) = ICx_side;
    params_used(L, 7) = ICx_width;
    params_used(L, 8) = FREQ_width;
    
    end
    
end
end

end
end
end

end
end
end

L = find(isnan(params_used(:, 1)), 1, 'first');
params_used(L:end, :) = [];

clear L

%% Limit the variables to plot if desired

params_used = params;

itds = 100;
repet = 1;

snr_iccl = 2;
snr_icx = 2;
icx_weight = 1;
icx_side = -0.5;
icx_width = -0.5;
freq_width = -0.5;

params_used = params_used(logical(sum(params_used(:, 1) == itds, 2)), :);
params_used = params_used(logical(sum(params_used(:, 2) == repet, 2)), :);
params_used = params_used(logical(sum(params_used(:, 3) == snr_iccl, 2)), :);
params_used = params_used(logical(sum(params_used(:, 4) == snr_icx, 2)), :);
params_used = params_used(logical(sum(params_used(:, 5) == icx_weight, 2)), :);
params_used = params_used(logical(sum(params_used(:, 6) == icx_side, 2)), :);
params_used = params_used(logical(sum(params_used(:, 7) == icx_width, 2)), :);
params_used = params_used(logical(sum(params_used(:, 8) == freq_width, 2)), :);

%% Or if you have a list of combinations which work
params_used = params;
ix = 108;
params_used = params_used(ix, :);

%% keep them all
params_used = params;
%%

% slopes = nan(size(params_used, 1), 3);
% intercepts = nan(size(params_used, 1), 3);

i = 61;

dir_name = 'K:\python\model_icx_codeonly\';

for p = 301:400%size(params_used, 1)
  
ITD = params_used(p, 1);
rep = params_used(p, 2);

SNR_ICcl = params_used(p, 3);
SNR_ICx = params_used(p, 4);
ICx_weight = params_used(p, 5);
    
ICx_side = params_used(p, 6);
ICx_width = params_used(p, 7);
FREQ_width = params_used(p, 8);

    
    r_txt = sprintf('%03d', rep);

    if isfile([dir_name 'SynchResults\data\testingvars_firstpass\modelout_ITD' num2str(ITD) '_' r_txt ...
    '_SNRICcl' num2str(SNR_ICcl * 10) ...
    '_SNRICx' num2str(SNR_ICx * 10) ...
    '_ICxweight' num2str(ICx_weight * 10) ...
    '_ICxside' num2str(ICx_side * 10) ...
    '_ICxwidth' num2str(ICx_width * 10) ...
    '_freqwidth' num2str(FREQ_width * 10) ...
    '_shiftshuff.mat'])

    disp('Parameters :')
    disp(['SNR_ICcl ' num2str(SNR_ICcl)])
    disp(['SNR_ICx ' num2str(SNR_ICx)])
    disp(['ICx_weight ' num2str(ICx_weight)])
    disp(['ICx_side ' num2str(ICx_side)])
    disp(['ICx_width ' num2str(ICx_width)])
    disp(['freq width ' num2str(FREQ_width)])
   
    eval(['load ' dir_name 'SynchResults\data\testingvars_firstpass\modelout_ITD' num2str(ITD) '_' r_txt ...
        '_SNRICcl' num2str(SNR_ICcl * 10) ...
        '_SNRICx' num2str(SNR_ICx * 10) ...
        '_ICxweight' num2str(ICx_weight * 10) ...
        '_ICxside' num2str(ICx_side * 10) ...
        '_ICxwidth' num2str(ICx_width * 10) ...
        '_freqwidth' num2str(FREQ_width * 10) ...
        '_shiftshuff.mat'])

figure

set(gcf, 'position', [0 0 700 800])

colors = {'k', 'r', 'b'};

subplot(4, 3, 4:5)
hold on

for t = 1:3

    b = syn(:, t);

    [slopes(p, t), intercepts(p, t)] = scatter_regression_plot(a(1:end), b(1:end), colors{t});
    xlabel('Geometric mean firing rate')
    ylabel({'Synchrony', '(coincidences/spike)'})

end
hold off


subplot(4, 3, 6)
binNum = 199;

sh = zeros(1, binNum);
ccg = zeros(1, binNum);
cor = zeros(1, binNum);

for ii = 1:size(synch_varying, 1)
    sh = nansum([sh; synch_varying{ii, 6}{1}']);
    ccg = nansum([ccg; synch_varying{ii, 5}{1}']);
    cor = nansum([cor; synch_varying{ii, 7}{1}']);
end

x = linspace(-100,100,binNum);

plot(x,ccg,'g');
hold on
plot(x,sh,'r');
plot(x,cor,'k');
hold off


subplot(4,3,1:3)
plot(bestITD,Rate,'o-')
ylabel('Firing rate (spikes/s)','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title(num2str(p))

subplot(4,3,7);
bsRaster(spikesICx(:,i),[0 T])
set(gca,'YTick',1:10:81,'YTickLabel',bestITD(1:10:81))
xlabel('Time (ms)','fontsize',15)
ylabel('Best ITD (\mus)','fontsize',15)
title('Original OT spikes','fontsize',15)

subplot(4,3,10);
bsRaster(spikes2(:,i),[0 T]);hold on
set(gca,'YTick',1:10:81,'YTickLabel',bestITD(1:10:81))
xlabel('Time (ms)','fontsize',15)
ylabel('Best ITD (\mus)','fontsize',15)
title('Modified OT spikes','fontsize',15)

% mean rate

Rm = mean(cellfun(@length,spikesICx(:,i)),2);
Rm2 = mean(cellfun(@length,spikes2(:,i)),2);

subplot(4,3,8);
plot(bestITD,Rm,'o-')
ylabel('Spike Count','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title('Original OT Spike Count','fontsize',15)

subplot(4,3,11);
plot(bestITD,Rm2,'o-')
ylabel('Spike Count','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title('Modified OT Spike Count','fontsize',15)


subplot(4,3,9);
plot(bestITD,SC,'o-')
ylabel('Spike Count','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title('Output Spike Count','fontsize',15)


subplot(4,3,12);
plot(bestITD,SC2,'o-')
ylabel('Spike Count','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title('Output Spike Count','fontsize',15)

% pause
    end
    
end

%%

y60 = 60 * slopes(:, 1) + intercepts(:, 1);
y0 = 0 * slopes(:, 1) + intercepts(:, 1);

x = [0 60];
OTvar = [1.3913e-05 8.5910e-05];

Y_OT = x * OTvar(1) + OTvar(2);
plot(x, Y_OT, 'k', 'LineWidth', 3);
hold on
for i = 1:length(y0)
    plot(x, [y0(i) y60(i)], 'r', 'LineWidth', 1)
end
plot([0 60], [9.2069e-04 9.2069e-04], 'b', 'LineWidth', 1, 'LineStyle', '--')
plot([0 60], [0 0], 'b', 'LineWidth', 1, 'LineStyle', '--')

xlabel('GM firing rate')
ylabel('synchrony')
title('Linear regressions for different parameters')
