function Synch4_get_ICx_Response_dynamicgain_Wweight_pipeline(repet, ITD, dir_name, SNR_ICcl, SNR_ICx, ICx_weight, ICx_side, ICx_width, freq_width, BW, gain_L, gain_U, W_weight)

%This file generates ICx spikes from ICcl spikes that are
%loaded in.
%It can take around 15 min.

tic

if ~exist('repet', 'var')
    repet = 0;
end

r_txt = sprintf('%03d', repet);

%% Simulation parameters

%Stimulus ITD                           %% Can change this value %%
if ~exist('ITD', 'var')
    ITD = 0;
end

%% Directory where SynchResults folder is found
if ~exist('dir_name', 'var')
    dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\'; %% Can change this value %%
end

eval(['cd ' dir_name '\SynchResults'])

%% Load parameters

eval(['load ' dir_name 'SynchResults\data\SynchParms_SNRICcl' num2str(SNR_ICcl * 10) ...
    '_SNRICx' num2str(SNR_ICx * 100000) ...
    '_ICxweight' num2str(ICx_weight * 10) ...
    '_ICxside' num2str(ICx_side * 10) ...
    '_ICxwidth' num2str(ICx_width * 10) ...
    '_freqwidth' num2str(freq_width * 10) ...
    '_bw' num2str(BW * 100000) ...
    '_Wweight' num2str(W_weight * 10)])

%Number of ICx neurons
NumN = size(We,1);

eval(['load ' dir_name 'SynchResults\data\cc\SynchCC_ITD100_001 Ts2'])


%% Load data

eval(['load ' dir_name 'SynchResults\data\iccl\SynchICcl_ITD' num2str(ITD) '_001_SNR' num2str(SNR_ICcl * 10) ' spikesICcl T BITD CF'])

[Nf,Nn,Ntrial] = size(spikesICcl);

%% Time

t = 0:Ts2:T;
Nt = length(t);

%% initialize

spikesICx = cell(NumN,Ntrial);

%%

for i = 1:Ntrial
  
    spikes = spikesICcl(:,:,i); 

    % Get ICx conductance
    [ge_icx,gi_icx] = getICxConductance(spikes,We,Wi,tau_syn,t); 
    
    
    %%%%%%%%% OLD NOISE CALCULATIONS
    % Added this reduction in input strength to attempt to reduce firing
    % rate range
%     
%     ge_icx = ICx_weight * ge_icx;   %% change
%     gi_icx =  gi_icx;   %% change
%     
%     %Get noise RMS values
%     rms_noise_e = sqrt(rms(ge_icx(:)))/SNR_ICx;
%         
%     rms_noise_i = sqrt(rms(gi_icx(:)))/SNR_ICx;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%% NEW noise calculations using a GAIN factor
    
    gain_e = rescale(mean(ge_icx, 2), gain_L, gain_U);

    rms_noise_e = (-1+sqrt(1+4.*SNR_ICx.^2.*gain_e))./(2.*SNR_ICx.^2);
    a_e = gain_e - rms_noise_e;
    ge_icx = ge_icx./rms(ge_icx(:)).*a_e;
    
    gain_i = rescale(mean(gi_icx, 2), gain_L, gain_U);

    rms_noise_i = (-1+sqrt(1+4.*SNR_ICx.^2.*gain_i))./(2.*SNR_ICx.^2);
    a_i = gain_i - rms_noise_i;
    gi_icx = gi_icx./rms(gi_icx(:)).*a_i;
    
        % Get ICx response
        for n = 1:NumN

            %Generate noise for E and I
            ne = genSignal(T+5,Ts2,rms_noise_e(n),2,2*pi*BW);
            ne = ne(1:Nt);
            
            ni = genSignal(T+5,Ts2,rms_noise_i(n),2,2*pi*BW);
            ni = ni(1:Nt);
        
            %Get spikes
            [~,fr,sp] = genSpikesICxAdExTrc(ge_icx(n,:) + ne,gi_icx(n,:) + ni,ge_icx(n,:)*0,Ts2,NP);
            
            sp = sp(1:fr);

            spikesICx{n,i} = sp;
                        
        end
        
end


%% Plots

% % plot rasters
% M = 20; %number of neurons to plot
% 
% figure(1);clf;
% for n = 1:M
%     subplot(M,1,n);
%     bsRaster(spikesICx(n*4+1,:).',[0 T])
%     
%     xlim([0 T])
%     if n < M
%         set(gca,'XTickLabels',{})
%     end
%     ylabel(num2str(bestITD(n*4+1)))
% end
% xlabel('Time (ms)','fontsize',15)
% 
% %Plot rate
Rate = squeeze(mean(cellfun(@length,spikesICx),2)*1000/T);
% 
% figure(2);clf;
% plot(bestITD,Rate,'o')
% xlabel('Best ITD (\mus)','fontsize',15)
% ylabel('Rate (spikes/s)','fontsize',15)


%% save data

eval(['save ' dir_name 'SynchResults\data\icx\SynchICx_ITD' num2str(ITD) '_' r_txt ...
    '_SNRICcl' num2str(SNR_ICcl * 10) ...
    '_SNRICx' num2str(SNR_ICx * 100000) ...
    '_ICxweight' num2str(ICx_weight * 10) ...
    '_ICxside' num2str(ICx_side * 10) ...
    '_ICxwidth' num2str(ICx_width * 10) ...
    '_freqwidth' num2str(freq_width * 10) ...
    '_bw' num2str(BW * 100000) ...
    '_gainL' num2str(gain_L * 100) ...
    '_gainU' num2str(gain_U) ...
    '_Wweight' num2str(W_weight * 10)])

toc

end
