FMIN_var = 0.25;
WEIGHT_var = 12;

dir_name = 'K:\python\model_icx_codeonly\';
% dir_name = 'F:\Desktop\Code\Reproducibility\model_ICx\';

for Fmin = FMIN_var
    f_txt = sprintf('%02d', Fmin * 10);
    for weight = WEIGHT_var

repet = 1;
load('K:\python\model_icx_codeonly\SynchResults\data\model_bounds.mat')
% load('f:\desktop\Code\Reproducibility\model_ICx\SynchResults\data\model_bounds.mat')
T = 150;
ix = find(itd == T);
itd = T;

n_rep = 10;

for t = 1:length(itd)
    tic
    spikesICx = cell(1, n_rep);
    spikes2 = cell(1, n_rep);
    output_orig = cell(1, n_rep);
    output_mani = cell(1, n_rep);
    for r = 1:n_rep
        disp(['running itd = ' num2str(itd(t)) '; rep = ' num2str(r)]) 
%         [spikesICx{r}, spikes2{r}, output_orig{1, r}, output_mani{1, r}] = decoder_calc(itd(t), repet, r, bounds(t, :), Fmin, weight, dir_name);
        [spikesICx{r}, spikes2{r}, output_orig{1, r}, output_mani{1, r}] = decoder_calc(itd(t), repet, r, bounds(ix, :), Fmin, weight, dir_name);
    end
    toc
    
    spikesICx = [spikesICx{1, :}];
    spikes2 = [spikes2{1, :}];
    bestITD = -200:5:200;
    
    save([dir_name 'SynchResults\data\tmp\SynchOut_ITD' num2str(itd(t)) '_' f_txt], 'spikesICx', 'spikes2', 'output_orig', 'output_mani', 'Fmin', 'weight', 'bounds')

end

    end
end
%%


function [spikesout1, spikesout2, SC, SC2] = decoder_calc(ITD, repet, trial, bounds, FMIN, weight, dir_name)
%%
if ~exist('repet', 'var')
    repet = 0;
end

r_txt = sprintf('%03d', repet);

% % Directory where SynchResults folder is found

if ~exist('dir_name', 'var')
%     dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\';
    dir_name = 'K:\python\model_icx_codeonly\';
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
    W(n,:) = weight*exp(-.5*((bestITD - bestITD(n))/30).^2); 
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

inter = 0:25:ceil(T);
Li = length(inter);

% % for each subinterval, shift spikes closer together.
g2 = g_out;

spikes2 = spikesICx;

%Get mean time in each subinterval
mean_time = zeros(1,Li-1);
for k = 1:Li-1
    
    sp_temp = [];
    for n = 1:L2
        sp2 = spikesICx{ind2(n),trial};
        sp_temp = [sp_temp sp2(sp2>=inter(k) & sp2<inter(k+1))];
    end
    mean_time(1, k) = mean(sp_temp);
end
%Shift times close to mean
for n = 1:L2
    sp2 = spikesICx{ind2(n),trial};
    
    sp_new = [];
    
    for k = 1:Li-1
        
        sp_int = sp2(sp2>=inter(k) & sp2<inter(k+1));
        
        d = sp_int - mean_time(1, k);
        
        %frac = rand; %Shift by a fraction of the distance. If frac = 1, they are all shifted to exactly the mean.
        fmin = FMIN;
        
        frac = fmin + (1 - fmin)*rand;
        
        sp_shift = sp_int - frac*d;
        
        sp_new = [sp_new sp_shift];
    end
    
    spikes2{ind2(n),trial} = sp_new;
    
    Nsp = length(sp_new);
    if Nsp>0
%         g2(ind2(n),:) = zeros(1, Nt);
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

V = zeros(Nn,Nt);

V2 = zeros(Nn,Nt);

for n = 1:Nn
    
    [V(n,:),SC(n)] = genSpikesICxAdExTrc(g_o(n,:),g(n,:)*0,g(n,:)*0,Ts2,NP);
    [V2(n,:),SC2(n)] = genSpikesICxAdExTrc(g(n,:),g(n,:)*0,g(n,:)*0,Ts2,NP);
    
end

spikesout1 = spikesICx(:, trial);
spikesout2 = spikes2(:, trial);

end
