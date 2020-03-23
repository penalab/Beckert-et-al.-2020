
% % Directory where SynchResults folder is found

if ~exist('dir_name', 'var')
    % dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\';  %% Can change this value %%
    dir_name = 'K:\python\model_icx_codeonly\';
end
eval(['cd ' dir_name '\SynchResults'])

% % Simulation parameters

%Stimulus ITD                           %% Can change this value %%
if ~exist('ITD', 'var')
    ITD = 80;
end
% % Load parameters

eval(['load ' dir_name 'SynchResults\data\iccl\SynchICcl_ITD' num2str(ITD) '_001 T'])
eval(['load ' dir_name 'SynchResults\data\SynchParms bestITD bestF NP'])
eval(['load ' dir_name 'SynchResults\data\cc\SynchCC_ITD' num2str(ITD) '_001 Ts2'])

% % Decoding population parameters %% Can change these values %%

%Synaptic time constant in ms
tau = 0.1;
%tau = 5;

% weights
Nn = length(bestITD);

W = zeros(Nn,Nn);
for n = 1:Nn
    %works with tau = .1
    W(n,:) = 12*exp(-.5*((bestITD - bestITD(n))/30).^2); 
    
    %works with tau = 5
    %W(n,:) = .25*exp(-.5*((bestITD - bestITD(n))/30).^2); %old
    
    %W(n,:) = .5*exp(-.5*((bestITD - bestITD(n))/30).^2);
end

% % Load ICx spikes

eval(['load ' dir_name 'SynchResults\data\icx\SynchICx_ITD' num2str(ITD) '_001 spikesICx Rate'])

% % plot mean rate

figure(1);clf;
plot(bestITD,Rate,'o-')
ylabel('Firing rate (spikes/s)','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)

% %

t = 0:Ts2:T;
Nt = length(t);

%% Get fitered spikes

Ntrial = size(spikesICx,2);

%Initialize filtered spike matrix
g_out = zeros(Nn,Nt);

%Pick a trial                               %% Can change this value %%
i = 2;

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

%% Find the side peak(s) and main peak %% Can change these values %%

%ITD = 150
%Find main peak
ind1 = find(bestITD >= 50);L1 = length(ind1);

%Find side peak
ind2 = find(bestITD <= -50);L2 = length(ind2);

%delay added to main peak spikes
T1 = 35;


%{
%Find main peak
ind1 = find(bestITD >= 55 & bestITD <= 190);L1 = length(ind1);

%Find side peak
ind2 = find(bestITD <= -55 & bestITD >= -190);L2 = length(ind2);

%delay added to main peak spikes
T1 = 0;
%}

% % time intevals

inter = 0:25:ceil(T);
Li = length(inter);


%% for each subinterval, shift spikes closer together.
g2 = g_out;

spikes2 = spikesICx;

%Get mean time in each subinterval
mean_time = zeros(1,Li-1);
for k = 1:Li-1
    
    sp_temp = [];
    for n = 1:L2
        sp2 = spikesICx{ind2(n),i};
        sp_temp = [sp_temp sp2(sp2>=inter(k) & sp2<inter(k+1))];
    end
    mean_time(1, k) = mean(sp_temp);
end
%Shift times close to mean
for n = 1:L2
    sp2 = spikesICx{ind2(n),i};
    
    sp_new = [];
    
    for k = 1:Li-1
        
        sp_int = sp2(sp2>=inter(k) & sp2<inter(k+1));
        
        d = sp_int - mean_time(1, k);
        
        %frac = rand; %Shift by a fraction of the distance. If frac = 1, they are all shifted to exactly the mean.
        fmin = 1;
        
        frac = fmin + (1 - fmin)*rand;
        
        sp_shift = sp_int - frac*d;
        
        sp_new = [sp_new sp_shift];
    end
    
    spikes2{ind2(n),i} = sp_new;
    
    Nsp = length(sp_new);
    if Nsp>0
        g2(ind2(n),:) = zeros(1, Nt);
        for l = 1:Nsp
            h = zeros(1,Nt);
            h(t >= sp_new(l)) = exp(-(t(t >= sp_new(l)) - sp_new(l))/tau);
            g2(ind2(n),:) = g2(ind2(n),:) + h;
        end
    end
    
end


% % Get output spikes

%Apply weights
g = W*g2;

g_o = W*g_out;

%Initialize spike counts
SC = zeros(1,Nn);

SC2 = zeros(1,Nn);

for n = 1:Nn
    
    [~,SC(n)] = genSpikesICxAdExTrc(g_o(n,:),g(n,:)*0,g(n,:)*0,Ts2,NP);
    [~,SC2(n)] = genSpikesICxAdExTrc(g(n,:),g(n,:)*0,g(n,:)*0,Ts2,NP);
    
end


%% Plot

figure(2);clf;
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
plot(bestITD,Rm,'o-')
ylabel('Spike Count','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title('Original OT Spike Count','fontsize',15)


subplot(2,3,5);
plot(bestITD,Rm2,'o-')
ylabel('Spike Count','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title('Modified OT Spike Count','fontsize',15)


subplot(2,3,3);
plot(bestITD,SC,'o-')
ylabel('Spike Count','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title('Output Spike Count','fontsize',15)


subplot(2,3,6);
plot(bestITD,SC2,'o-')
ylabel('Spike Count','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)
title('Output Spike Count','fontsize',15)


