function Synch3_get_ICcl_Response(repet, ITD, dir_name)

tic

if ~exist('repet', 'var')
    repet = 0;
end

r_txt = sprintf('%03d', repet);

%This file generates ICcl spikes from cross correlation input that is
%loaded in.
%It can take around 15 min.

%% Simulation parameters

%Stimulus ITD                           %% Can change this value %%
if ~exist('ITD', 'var')
    ITD = 0;
end

%Number of trials                       %% Can change this value %%
Ntrial = 100;

%% Directory where SynchResults folder is found
if ~exist('dir_name', 'var')
    dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\'; %% Can change this value %%
end

eval(['cd ' dir_name '\SynchResults'])

%% Load parameters

eval(['load ' dir_name 'SynchResults/data/SynchParms_reducedSNR s BW SNR_ICcl BITD tau_syn NP'])

%% Load cross correlation inputs at given ITD

eval(['load ' dir_name 'SynchResults/data/cc/SynchCC_ITD' num2str(ITD) '_' r_txt ' CF Ts Ts2 X'])

%Size of input
[Nf,Nd,Nt] = size(X);

%% Stimulus time

T = (Nt-1)*Ts2;
t = 0:Ts2:T;

%% Get ICcl spikes

[spikesICcl,R] = getICclAdExSpikesTrials(X,BITD,Ntrial,Ts2,T,SNR_ICcl,BW,s,NP);

%% save

eval(['save ' dir_name 'SynchResults/data/iccl/SynchICcl_ITD' num2str(ITD) '_' r_txt '_reducedSNR R spikesICcl SNR_ICcl BW s NP T BITD CF'])


%% Plots

%Raster plots

% figure(1);clf;
% for n = 1:9
%     subplot(3,3,n);
%     bsRaster(squeeze(spikesICcl(n+30,33,:)),[0 T])
%     xlim([0 T])
%     if n == 2
%     title('Example ICcl Rasters')
%     end
%     xlabel('Time (ms)','fontsize',15)
%     ylabel('Trials','fontsize',15)
% 
% end
% 
% % Spike count
% Rm = mean(R,3);
% 
% figure(2);clf;
% imagesc(BITD,CF,Rm);axis xy;colorbar
% xlabel('Best ITD (\mus)','fontsize',15)
% ylabel('Best Frequency (Hz)','fontsize',15)
% title('ICcl Spike Counts','fontsize',15)

toc
