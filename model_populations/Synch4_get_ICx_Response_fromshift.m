function Synch4_get_ICx_Response_fromshift(max_shift, dir_name)

%This file generates ICx spikes from ICcl spikes that are
%loaded in.
%It can take around 15 min.

tic

if ~exist('max_shift', 'var')
    max_shift = 100;
end

%% Simulation parameters

%Stimulus ITD                           %% Can change this value %%
if ~exist('ITD', 'var')
    ITD = 100;
end

%% Directory where SynchResults folder is found
if ~exist('dir_name', 'var')
    dir_name = 'K:\python\model_icx_codeonly\'; %% Can change this value %%
end

eval(['cd ' dir_name '\SynchResults'])

%% Load parameters

eval(['load ' dir_name 'SynchResults\data\SynchParms bestITD bestF We Wi tau_syn BW NP SNR_ICx '])

%Number of ICx neurons
NumN = size(We,1);

eval(['load ' dir_name 'SynchResults\data\cc\SynchCC_ITD' num2str(100) '_001 Ts2'])


%% Load data

eval(['load ' dir_name 'SynchResults\model_synch_iccl_shifttest_' num2str(max_shift) ' spikesICcl_shift T BITD CF'])

spikesICcl = spikesICcl_shift;
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
    
    %Get noise RMS values
    rms_noise_e = sqrt(rms(ge_icx(:)))/SNR_ICx;
        
    rms_noise_i = sqrt(rms(gi_icx(:)))/SNR_ICx;
        
        % Get ICx response
        for n = 1:NumN

            %Generate noise for E and I
            ne = genSignal(T+5,Ts2,rms_noise_e,2,2*pi*BW);
            ne = ne(1:Nt);
            
            ni = genSignal(T+5,Ts2,rms_noise_i,2,2*pi*BW);
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

eval(['save ' dir_name 'SynchResults\model_synch_icx_shifttest_' num2str(max_shift) ' spikesICx Rate'])

toc

end
