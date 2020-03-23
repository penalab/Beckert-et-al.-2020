%% Run once per parameter sets

cd('K:\python\model_icx_codeonly\SynchResults')
Synch1_Model_Parameters_ReducedWeights

% % Iterate across variables

% itds = -200:20:200;
itds = 100;
repet = 1;

for rep = repet
    for ITD = itds
%         disp(['calculating cross-correlation : rep = ' num2str(rep) '; ITD = ' num2str(ITD)])
%         Synch2_get_CC_Response(rep, ITD, 'K:\python\model_icx_codeonly\')
%         disp(['calculating ICcl responses : rep = ' num2str(rep) '; ITD = ' num2str(ITD)])
%         Synch3_get_ICcl_Response_reducedSNR(rep, ITD, 'K:\python\model_icx_codeonly\', 2)
        disp(['calculating ICx responses : rep = ' num2str(rep) '; ITD = ' num2str(ITD)])
        Synch4_get_ICx_Response_reducedSNR(rep, ITD, 'K:\python\model_icx_codeonly\', 2)
        Synch4b_adjust_ICx_Response_fromshift_RunNCalc
        Synch5_Decode_ICxvShiftTimes_trials
    end
end

%% Plotting

dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\'; %% Can change this value %%

ITD = -100;
rep = 1;
r_txt = sprintf('%03d', rep);

load([dir_name 'SynchResults\data\SynchParms'])
load([dir_name 'SynchResults\data\cc\SynchCC_ITD' num2str(ITD) '_' r_txt])
load([dir_name 'SynchResults\data\iccl\SynchICcl_ITD' num2str(ITD) '_' r_txt])

figure(1);clf;
for n = 1:9
    subplot(3,3,n);
    bsRaster(squeeze(spikesICcl(n+30,33,:)),[0 T]);
    xlim([0 T])
    if n == 2
    title('Example ICcl Rasters')
    end
    xlabel('Time (ms)','fontsize',15)
    ylabel('Trials','fontsize',15)
end

% Spike count
Rm = mean(R,3);

figure(2);clf;
imagesc(BITD,CF,Rm);axis xy;colorbar
xlabel('Best ITD (\mus)','fontsize',15)
ylabel('Best Frequency (Hz)','fontsize',15)
title('ICcl Spike Counts','fontsize',15)

load([dir_name 'SynchResults\data\icx\SynchICx_ITD' num2str(ITD) '_' r_txt])

% plot rasters
M = 20; %number of neurons to plot

figure(3);clf;
for n = 1:M
    subplot(M,1,n);
    bsRaster(spikesICx(n*4+1,:).',[0 T])
    
    xlim([0 T])
    if n < M
        set(gca,'XTickLabels',{})
    end
    ylabel(num2str(bestITD(n*4+1)))
end
xlabel('Time (ms)','fontsize',15)

%Plot rate
Rate = squeeze(mean(cellfun(@length,spikesICx),2)*1000/T);

figure(4);clf;
plot(bestITD,Rate,'o')
xlabel('Best ITD (\mus)','fontsize',15)
ylabel('Rate (spikes/s)','fontsize',15)

