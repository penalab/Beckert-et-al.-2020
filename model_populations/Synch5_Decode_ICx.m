repet = 1;
ITD = 100;

tic

if ~exist('repet', 'var')
    repet = 0;
end

r_txt = sprintf('%03d', repet);

% % Directory where SynchResults folder is found

% dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\'; %% Can change this value %%
dir_name = 'K:\python\model_icx_codeonly\';

eval(['cd ' dir_name '\SynchResults'])

% % Simulation parameters

%Stimulus ITD                           %% Can change this value %%
if ~exist('ITD', 'var')
    ITD = 0;
end

% % Load parameters

eval(['load ' dir_name 'SynchResults\data\iccl\SynchICcl_ITD' num2str(ITD) '_' r_txt ' T'])
eval(['load ' dir_name 'SynchResults\data\SynchParms bestITD bestF NP'])
% eval(['load ' dir_name 'SynchResults\data\cc\SynchCC_ITD' num2str(ITD) '_' r_txt ' Ts2'])
Ts2 = 0.0205;

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
    %W(n,:) = 3.5*exp(-.5*((bestITD - bestITD(n))/40).^2); 
    
    %works with tau = 5
    %W(n,:) = .25*exp(-.5*((bestITD - bestITD(n))/30).^2);
end

% % Load ICx spikes

eval(['load ' dir_name 'SynchResults\data\icx\SynchICx_ITD' num2str(ITD) '_' r_txt ' spikesICx Rate'])

% % plot mean rate

figure(1);clf;
plot(bestITD,Rate,'o-')
ylabel('Firing rate (spikes/s)','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)

% %

t = 0:Ts2:T;
Nt = length(t);

% % Get fitered spikes

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


% %Find main peak
% ind1 = find(bestITD >= -190 & bestITD <= -75);
% L1 = length(ind1);
% 
% %Find side peak
% ind2 = find(bestITD <= 185 & bestITD >= 85);
% L2 = length(ind2);
% 
% %delay added to main peak spikes
% T1 = 0;


%% Substitute main peak spikes into side peak(s)

g2 = g_out;

spikes2 = spikesICx;

% for r = 1:Ntrial
for n = 1:L2
    
    g2(ind2(n),:) = g2(ind2(n),:)*0;
    % get matching neurons
    
    n1 = 0;
    sp2 = spikesICx{ind2(n),i};
    
    nsp = length(sp2);
    
    if nsp>0

        cnt = n;
        
        while n1 < nsp

        sp1 = spikesICx{ind1(cnt),i};
        n1 = length(sp1);
        if n < L1/2
        cnt = cnt + 1;
        else
        cnt = cnt - 1;
        end

        end
      

        sp = sp1(1:nsp) + T1;

        spikes2{ind2(n),i} = sp;

        Nsp = length(sp);
        if Nsp>0
            for l = 1:Nsp
                h = zeros(1,Nt);
                h(t >= sp(l)) = exp(-(t(t >= sp(l)) - sp(l))/tau);
                g2(ind2(n),:) = g2(ind2(n),:) + h;
            end
        end
    
    end
end
% end

% % Get output spikes

%Apply weights
g = W*g2;

g_o = W*g_out;

%Initialize spike counts
SC = zeros(1,Nn);

SC2 = zeros(1,Nn);

V = zeros(Nn,Nt);

V2 = zeros(Nn,Nt);

for n = 1:Nn
    
    [V(n,:),SC(n)] = genSpikesICxAdExTrc(g_o(n,:),g(n,:)*0,g(n,:)*0,Ts2,NP);
    [V2(n,:),SC2(n)] = genSpikesICxAdExTrc(g(n,:),g(n,:)*0,g(n,:)*0,Ts2,NP);
    
end

%%

eval(['save ' dir_name 'SynchResults\data\out\SynchOut_ITD' num2str(ITD) '_' r_txt ' spikesICx spikes2 SC SC2'])


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

%%
figure(3);clf;
subplot(211)
imagesc(t,bestITD,V);axis xy;colorbar

subplot(212)
imagesc(t,bestITD,V2);axis xy;colorbar

%% This is a test plot to confirm that the spike trains were shifted correctly

spk = spikesICx;
spk2 = spikes2;

for i = 1:size(spk, 1)

    figure
    y = get_spike_y(spk(i, :));
    scatter(cell2mat(spk(i, :)), cell2mat(y), 10, 'fill')
    hold on
    y = get_spike_y(spk2(i, :));
    scatter(cell2mat(spk2(i, :)), cell2mat(y), 30)
    title(num2str(i))
    
    pause
    close all
end