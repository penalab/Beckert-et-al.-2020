function Synch2_get_CC_Response(repet, ITD, dir_name)

if ~exist('repet', 'var')
    repet = 0;
end

r_txt = sprintf('%03d', repet);

tic

%This runs the cross-correlation model for one ITD. It can take around 5
%minutes.

%% Directory where SynchResults folder is found

if ~exist('dir_name', 'var')
    dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\'; %% Can change this value %%
end

eval(['cd ' dir_name '\SynchResults'])

%% Stimulus parameters
%ITD,ILD,ABL %% Can change these values %%

if ~exist('ITD', 'var')
    ITD=0;
end

ILD=0;
ABL=50;

%Noise (standard deviation = sig*(current signal value))
sig=0;

%% Time

%Time step for original sampling rate (ms)      %% Can change this value %%
Ts2 = 1000/48828; 

%Upsampling factor before cross correlation     %% Can change this value %%
ds = 20;

%New time step (ms)
Ts = Ts2/ds;

%% Load frequencies

eval(['load ' dir_name 'SynchResults/data/SynchParms CF'])

Nf = length(CF);

%% rms values

RMSL=ABL-.5*ILD;
RMSR=ABL+.5*ILD;
rmsL=2*10^(RMSL/20);
rmsR=2*10^(RMSR/20);

%% Load sound from figure 3

%Get sound
load('006-2015-0302-02-FrozenILD.mat')
s = curvesettings.stimcache.S0;

%set RMS to impose ILD and ABL
sL = rmsL*s(1,:)/rms(s(1,:));
sR = rmsR*s(1,:)/rms(s(2,:));

%Upsample
NL = resample(sL,ds,1);
NR = resample(sR,ds,1);

%% Add ITD

%Get a buffer for shifting sounds
Lt=length(NL);
ind=300:Lt-299;

%indices of shifted sound
in=ITD/1000/Ts;

%Shift to impose ITD
sL=NL(round(ind+in));
sR=NR(ind);

%% Generate Gammatone filter coefficients

FS=1/Ts*1000;
fcoefs=getFilterCoefs(CF,FS);

%% Run simulation

%%%%%%%%%%%%%%Filter Bank%%%%%%%%%%%%%%%%%

vL=ERBFilterBank(sL,fcoefs);
vR=ERBFilterBank(sR,fcoefs);

%Noise
if sig > 0
vL=vL+abs(vL).*randn(size(vL))*sig;
vR=vR+abs(vR).*randn(size(vR))*sig;
end

clear sL sR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%Cross correlation%%%%%%%%%%%%%%
%cross correlation in first frequency channel
x=CrossCorrelationwGainControl(vL(1,:),vR(1,:),Ts,sig);

%Get internal delays
ITDin=linspace(-200,200,size(x,2));

%Initialize
X = zeros(Nf,length(ITDin),size(x(:,:,1:ds:end),3));

%First frequency channel. Downsample
X(1,:,:) = x(:,:,1:ds:end);

%Run for the other frequencies
for n = 2:Nf
    
    x=CrossCorrelationwGainControl(vL(n,:),vR(n,:),Ts,sig);
    
    X(n,:,:) = x(:,:,1:ds:end); %Downsample cross correlation
end

clear vL vR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Save

eval(['save ' dir_name 'SynchResults/data/cc/SynchCC_ITD' num2str(ITD) '_' r_txt ' CF Ts Ts2 X ITD ILD ABL'])

toc

end