clear;

%% Directory where SynchResults folder is found

% dir_name = '/Users/fischer9/Dropbox/MATLAB/ReliabilitySynchMike/';  %% Can change this value %%
% 
% eval(['cd ' dir_name '/SynchResults'])

%% Simulation parameters

fmin = .25;

%Stimulus ITD                           %% Can change this value %%
%ITD = 150;

%% Load parameters

eval(['load K:\python\model_icx_codeonly\SynchResults\data\testingoutputfrombriansScript\SynchICcl_ITD' num2str(150) ' T'])
load('K:\python\model_icx_codeonly\SynchResults\data\testingoutputfrombriansScript\SynchParms', 'bestITD', 'bestF', 'NP')
eval(['load K:\python\model_icx_codeonly\SynchResults\data\testingoutputfrombriansScript\SynchCC_ITD' num2str(150) ' Ts2'])

load('K:\python\model_icx_codeonly\SynchResults\data\model_bounds.mat')

%Litd = length(itd);

IND = [2:7 15:20];
Litd = length(IND);

%% Fix bounds

bounds(IND(1),[3 4]) = [57 79];
bounds(IND(2),[3 4]) = [55 79];
bounds(IND(3),4) = 79;
bounds(IND(4),[3 4]) = [56 79];
bounds(IND(5),[3 4]) = [48 79];
bounds(IND(6),[3 4]) = [56 78];
bounds(IND(7),4) = 25;
bounds(IND(8),[3 4]) = [3 31];
bounds(IND(9),[3 4]) = [2 26];
bounds(IND(10),[3 4]) = [3 25];
bounds(IND(11),[3 4]) = [3 29];
bounds(IND(12),[3 4]) = [2 30];

%%

t = 0:Ts2:T;
Nt = length(t);

%% time intevals

inter = 0:25:ceil(T);
Li = length(inter);

%% Decoding population parameters %% Can change these values %%

%Synaptic time constant in ms
tau = 0.1;
%tau = 5;

% weights
Nn = length(bestITD);

W = zeros(Nn,Nn);
for n = 1:Nn
    %works with tau = .1
    W(n,:) = 11*exp(-.5*((bestITD - bestITD(n))/30).^2); %works
    
    %W(n,:) = 25*exp(-.5*((bestITD - bestITD(n))/10).^2); %works
    %W(n,:) = 50*exp(-.5*((bestITD - bestITD(n))/5).^2); %works
    %W(n,:) = 7*exp(-.5*((bestITD - bestITD(n))/75).^2); %works
    
    
    
    %works with tau = 5
    %W(n,:) = .25*exp(-.5*((bestITD - bestITD(n))/30).^2); %old
    
    %W(n,:) = .5*exp(-.5*((bestITD - bestITD(n))/30).^2);
end

%%
Ntrial = 100;

%Initialize spike counts
SC = zeros(Litd,Ntrial,Nn);

SC2 = zeros(Litd,Ntrial,Nn);

%%

for m = 8

%% Load ICx spikes

%eval(['load /Users/fischer9/Dropbox/icx/SynchICx_ITD' num2str(itd(m)) '_001.mat'])
%eval(['load /Users/fischer9/Dropbox/icx/SynchICx_ITD' num2str(itd(IND(m))) '_001.mat'])

% load('K:\python\model_icx_codeonly\SynchResults\data\testingoutputfrombriansScript\SynchICx_ITD100_001_SNRICcl20_SNRICx2000000_ICxweight5_ICxside-5_ICxwidth-5_freqwidth-5_bw10000_gainL100_gainU1_shiftshuff.mat','spikesICx')
load('K:\python\model_icx_codeonly\SynchResults\data\icx\SynchICx_ITD100_001_SNRICcl20_SNRICx1000000_ICxweight20_ICxside-5_ICxwidth-5_freqwidth-5_bw10000_gainL100_gainU1_shiftshuff.mat', 'spikesICx')
%% Find the side peak(s) %% Can change these values %%

ind2 = find(bestITD <= bestITD(bounds(IND(m),4)) & bestITD >= bestITD(bounds(IND(m),3)));L2 = length(ind2);


%% Get fitered spikes

spikes2 = spikesICx;

for i = 1:Ntrial

%%


%Initialize filtered spike matrix
g_out = zeros(Nn,Nt);

%Pick a trial                               %% Can change this value %%
%i = 2;

for n = 1:Nn
    %Get spikes for this neuron
    sp = spikesICx{n,i};
        
    %Add exponential at each time
    Nsp = length(sp);
    if Nsp>0
        for l = 1:Nsp
            h = zeros(1,Nt);
            h(t >= sp(l)) = exp(-(t(t >= sp(l)) - sp(l))/tau);
            g_out(n,:) = g_out(n,:) + h;
        end
    end
end


%% for each subinterval, shift spikes closer together.
g2 = g_out;

%spikes2 = spikesICx;

%Get mean time in each subinterval
mean_time = zeros(1,Li-1);

for k = 1:Li-1
    
    sp_temp = [];
    for n = 1:L2
        sp2 = spikesICx{ind2(n),i};
        sp_temp = [sp_temp sp2(sp2>=inter(k) & sp2<inter(k+1))];
    end
    mean_time(k) = mean(sp_temp);
end
    
%Shift times close to mean
for n = 1:L2
    sp2 = spikesICx{ind2(n),i};
    
    sp_new = [];
    
    for k = 1:Li-1
        
        sp_int = sp2(sp2>=inter(k) & sp2<inter(k+1));
        
        d = sp_int - mean_time(k);
        
        %frac = rand*.5; %Shift by a fraction of the distance. If frac = 1, they are all shifted to exactly the mean.
                
        frac = fmin + (1 - fmin)*rand;
        
        sp_shift = sp_int - frac*d;
        
        sp_new = [sp_new sp_shift];
    end
    
    spikes2{ind2(n),i} = sp_new;
    
    Nsp = length(sp_new);
        if Nsp>0
            for l = 1:Nsp
                h = zeros(1,Nt);
                h(t >= sp_new(l)) = exp(-(t(t >= sp_new(l)) - sp_new(l))/tau);
                g2(ind2(n),:) = g2(ind2(n),:) + h;
            end
        end
    
end

%% Get output spikes

%NP(13) = -60;NP(16) = -70;
NP(13) = -57.5;NP(16) = -67.5;

%Apply weights
g = W*g2;

g_o = W*g_out;


for n = 1:Nn
    
    %[~,SC(n)] = genSpikesICxAdExTrc(g_o(n,:),g(n,:)*0,g(n,:)*0,Ts2,NP);
    %[~,SC2(n)] = genSpikesICxAdExTrc(g(n,:),g(n,:)*0,g(n,:)*0,Ts2,NP);

    [~,SC(m,i,n)] = genSpikesICxAdExTrc(g_o(n,:),g(n,:)*0,g(n,:)*0,Ts2,NP);
    [~,SC2(m,i,n)] = genSpikesICxAdExTrc(g(n,:),g(n,:)*0,g(n,:)*0,Ts2,NP);

end


end


%% Plot

figure(IND(m));clf;
subplot(2,3,1);
bsRaster(spikesICx(:,i),[0 T])
set(gca,'YTick',1:10:81,'YTickLabel',bestITD(1:10:81))
xlabel('Time (ms)','fontsize',15)
ylabel('Best ITD (\mus)','fontsize',15)
title('Original OT spikes','fontsize',15)

subplot(2,3,4);
bsRaster(spikes2(:,i),[0 T]);hold on
set(gca,'YTick',1:10:81,'YTickLabel',bestITD(1:10:81))
xlabel('Time (ms)','fontsize',15)
ylabel('Best ITD (\mus)','fontsize',15)
title('Modified OT spikes','fontsize',15)

% mean rate

Rm = mean(cellfun(@length,spikesICx(:,i)),2);
Rm2 = mean(cellfun(@length,spikes2(:,i)),2);



subplot(2,3,2);
plot(bestITD,Rm,'o-');hold on
plot(bestITD(ind2),Rm(ind2),'ro')
ylabel('Spike Count','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title('Original OT Spike Count','fontsize',15)


subplot(2,3,5);
plot(bestITD,Rm2,'o-')
ylabel('Spike Count','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title('Modified OT Spike Count','fontsize',15)


M = squeeze(mean(SC,2));
M2 = squeeze(mean(SC2,2));

subplot(2,3,3);
%plot(bestITD,SC,'o-')
plot(bestITD,M,'o-')

ylabel('Spike Count','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title('Output Spike Count','fontsize',15)



subplot(2,3,6);
%plot(bestITD,SC2,'o-')


plot(bestITD,M2,'o-')

ylabel('Spike Count','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title('Output Spike Count','fontsize',15)


figure(100);clf;
plot(bestITD,Rm,'o-');hold on
plot(bestITD(ind2),Rm(ind2),'ro')
ylabel('Spike Count','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title('Original OT Spike Count','fontsize',15)



%SPH(m) = max(SC2(ind2))/max(SC2);

%%
clear sc sc2

sc(:,:) = SC(m,:,:);
sc2(:,:) = SC2(m,:,:);

%eval(['save ' dir_name 'SynchResults/ManipulateICx_ITD' num2str(itd(IND(m))) ' spikes2 sc sc2'])





end


%% plot
M = mean(SC,2);
M2 = mean(SC2,2);
    
SPH = zeros(1,Litd);

figure(1000);clf;
%for m = 1:Litd
    
    ind2 = find(bestITD <= bestITD(bounds(IND(m),4)) & bestITD >= bestITD(bounds(IND(m),3)));L2 = length(ind2);

    SPH(m) = max(M2(m,ind2))/max(M2(m,:));
    %subplot(3,4,m)
    plot(bestITD,M(m,:),'ro-',bestITD,M2(m,:),'bs-');
    ylabel('Spike Count','fontsize',15)
    xlabel('Best ITD (\mus)','fontsize',15)
    
    title(['ITD = ' num2str(itd(IND(m))) ', SPH = ' num2str(SPH(m))])
%end
   
    

%%



%%

%eval(['save ' dir_name 'SynchResults/SynchDecode bounds itd bestITD W tau SC SC2 IND fmin'])












