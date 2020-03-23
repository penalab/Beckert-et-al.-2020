clear;

%% Directory where SynchResults folder is found

% dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\';
dir_name = 'K:\python\model_icx_codeonly\';

%% CFs

CF = 1000:200:9000;

Nf = length(CF);

%% ICcl parameters

%Gain
s = 45;

%Noise bandwidth (kHz)
BW = 1;

%Signal to noise ratio
SNR_ICcl = 3;

%Number of ICcl neurons per frequency
Nn = 65;

%Best ITD of ICcl neurons (us)
BITD = (-(Nn-1)/2:(Nn-1)/2)*6;

%% ICx population parameters

%Best ITD of ICx
bestITD = -200:5:200; 

%Number of ICx neurons
NumN = length(bestITD);

%ITD weight SD of ICx neurons (us)
sITD = 10;

%Best frequencies of ICx neurons (Hz)
bestF = 6500 - abs(bestITD)/200*4000;

%Frequency weight SD of ICx neurons (Hz)
sF = 2500 - abs(bestF)/9000*1000;

%Signal to noise ratio
SNR_ICx = 2;

%% ICcl - ICx connection & synapse parameters

%Half-period
half_p = 0.5./bestF*1e6;

%Weight matrix
W = zeros(NumN,Nf,Nn);

%Weights are Gabor-like
for k = 1:NumN

    for n = 1:Nf
        %Frequency weights
        Wf = exp(-.5*((bestF(k) - CF(n))/sF(k)).^2);
        %ITD weights
        Witd = exp(-.5*((bestITD(k) - BITD)/sITD).^2) - .5*exp(-.5*((bestITD(k)+half_p(k) - BITD)/(sITD)).^2) - .5*exp(-.5*((bestITD(k)-half_p(k) - BITD)/(sITD)).^2);
        %Product
        W(k,n,:) = Wf*Witd;
    end

end

%Excitatory weights
We = W.*(W>0)*.032;

%Inhibitory weights
Wi = -W.*(W<0)*.024;

%Synaptic time constant (ms)
tau_syn = 5;

%% Biophysical parameters

Tref = 1;
Vrest = -70;
Vreset = -55;
EL = -75;
tRC = 60;
C = 5;
Ee = 0;
Ei = -80;
th_range = 0;
Del = 1;

%Adaptive threshold parameters
tau_th = 5;
alpha = 0;
ka = 5;
VT = -50;
Vi = -60;
ki = 1;

NP = [Tref, Vrest, Vreset, EL, tRC, C, Ee, Ei, th_range, Del, tau_th, alpha, VT, ka, ki, Vi];

%% Save

eval(['save ' dir_name 'SynchResults/data/SynchParms CF BITD s BW SNR_ICcl bestITD bestF We Wi SNR_ICx NP tau_syn'])
