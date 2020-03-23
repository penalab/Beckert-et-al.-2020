
repet = 1;
% load('K:\python\model_icx_codeonly\SynchResults\data\model_bounds.mat')
load('f:\desktop\Code\Reproducibility\model_ICx\SynchResults\data\model_bounds.mat')

n_rep = 100;

output_orig = cell(length(itd), n_rep);
output_mani = cell(length(itd), n_rep);

for t = 1:length(itd)
    tic
    for r = 1:n_rep
        disp(['running itd = ' num2str(itd(t)) '; rep = ' num2str(r)]) 
        [output_orig{t, r}, output_mani{t, r}] = decoder_calc(itd(t), repet, r, bounds(t, :));
    end
    toc
end

%%


function [SC, SC2] = decoder_calc(ITD, repet, trial, bounds, dir_name)
%%
if ~exist('repet', 'var')
    repet = 0;
end

r_txt = sprintf('%03d', repet);

% % Directory where SynchResults folder is found

if ~exist('dir_name', 'var')
    dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\';
    % dir_name = 'K:\python\model_icx_codeonly\';
end

eval(['cd ' dir_name '\SynchResults'])

% % Simulation parameters

%Stimulus ITD                           %% Can change this value %%
if ~exist('ITD', 'var')
    ITD = 0;
end

% % Load parameters

eval(['load ' dir_name 'SynchResults\data\iccl\SynchICcl_ITD' num2str(ITD) '_' r_txt ' T'])
eval(['load ' dir_name 'SynchResults\data\SynchParms bestITD NP'])
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
load([dir_name 'SynchResults\data\icx\SynchICx_ITD' num2str(ITD) '_' r_txt], 'spikesICx')

% %

t = 0:Ts2:T;
Nt = length(t);

%Initialize filtered spike matrix
g_out = zeros(Nn,Nt);


for n = 1:Nn
    %Get spikes for this neuron
    sp = spikesICx{n,trial};
        
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


%Find main peak
ind1 = find(bestITD >= bestITD(bounds(1)) & bestITD <= bestITD(bounds(2)));L1 = length(ind1);

%Find side peak
ind2 = find(bestITD <= bestITD(bounds(4)) & bestITD >= bestITD(bounds(3)));L2 = length(ind2);

%delay added to main peak spikes
T1 = 0;

% % Substitute main peak spikes into side peak(s)

g2 = g_out;

spikes2 = spikesICx;

for n = 1:L2
    
    g2(ind2(n),:) = g2(ind2(n),:)*0;
    % get matching neurons
    
    n1 = 0;
    sp2 = spikesICx{ind2(n),trial};
    
    nsp = length(sp2);
    
    if nsp>0

        cnt = n;
        
        while n1 < nsp

        sp1 = spikesICx{ind1(cnt),trial};
        n1 = length(sp1);
        if n < L1/2
        cnt = cnt + 1;
        else
        cnt = cnt - 1;
        end

        end
      

        sp = sp1(1:nsp) + T1;

        spikes2{ind2(n),trial} = sp;

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

end
