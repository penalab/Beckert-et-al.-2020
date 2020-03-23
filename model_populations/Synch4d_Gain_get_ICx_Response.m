clear;

%This file generates ICx spikes from ICcl spikes that are
%loaded in.
%It can take around 15 min.

%% Simulation parameters

%Stimulus ITD                           %% Can change this value %%
ITD = 100;

%% Directory where SynchResults folder is found

dir_name = '/Users/fischer9/Dropbox/MATLAB/ReliabilitySynchMike/';  %% Can change this value %%

eval(['cd ' dir_name '/SynchResults'])

%% Load parameters

eval(['load ' dir_name 'SynchResults/SynchParms bestITD bestF We Wi tau_syn BW NP SNR_ICx '])

%Number of ICx neurons
NumN = size(We,1);

eval(['load ' dir_name 'SynchResults/SynchCC_ITD' num2str(ITD) ' Ts2'])


%% Load data

eval(['load ' dir_name 'SynchResults/SynchICcl_ITD' num2str(ITD) ' spikesICcl T BITD CF'])

[Nf,Nn,Ntrial] = size(spikesICcl);

%% Time

t = 0:Ts2:T;
Nt = length(t);

%% initialize

%Ntrial = 10;

spikesICx = cell(NumN,Ntrial);

%% set the sum of sqrt(rms(ge_icx(:))) and rms_noise_e

gain = 10;
SNR_ICx = .2;

%%

for i = 1:Ntrial
  
    spikes = spikesICcl(:,:,i); 

    % Get ICx conductance
    [ge_icx,gi_icx] = getICxConductance(spikes,We,Wi,tau_syn,t); 
    
    %Get noise RMS values
    %rms_noise_e = sqrt(rms(ge_icx(:)))/SNR_ICx;
        
    %rms_noise_i = sqrt(rms(gi_icx(:)))/SNR_ICx;
    
    rms_noise_e = (-1+sqrt(1+4*SNR_ICx^2*gain))/(2*SNR_ICx^2);
    a_e = gain - rms_noise_e;
    ge_icx = ge_icx/rms(ge_icx(:))*a_e;
    
    rms_noise_i = (-1+sqrt(1+4*SNR_ICx^2*gain))/(2*SNR_ICx^2);
    a_i = gain - rms_noise_i;
    gi_icx = gi_icx/rms(gi_icx(:))*a_i;
    
        
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

% plot rasters
M = 20; %number of neurons to plot

figure(1);clf;
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

figure(2);clf;
plot(bestITD,Rate,'o')
xlabel('Best ITD (\mus)','fontsize',15)
ylabel('Rate (spikes/s)','fontsize',15)


%% save data

%eval(['save ' dir_name 'SynchResults/SynchICx_ITD' num2str(ITD) ' spikesICx Rate'])


