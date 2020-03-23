function [sL,sR]=genSignal_ITDILD(ITD,ILD,ABL,T,Ts)


%[sL,sR]=genSignal_ITDILD(ITD,ILD,ABL,T,Ts);
%This function generates noise inputs with the desired
%ITD, ILD, and ABL.
%ITD in us, must be multiple of time step
%ILD in dB
%ABL in dB
%T in ms
%Ts in ms
%Assumes a noise signal with 0 - 12 kHz, on and off ramps of 5 msec


t=0:Ts:T;
Lt=length(t);
%ind=150:Lt+149;
ind=300:Lt+299;
%ind=2000:Lt+1999;

[j,rampInd]=min(abs(t-5));

bandwidth=2*pi*12;

RMSL=ABL-.5*ILD;
RMSR=ABL+.5*ILD;
rmsL=2*10^(RMSL/20);
rmsR=2*10^(RMSR/20);

s=genSignal(T+10,Ts,1,2,bandwidth);

in=ITD/1000/Ts;
if in~=round(in)
    disp('Error: ITD not multiple of time step')
end
  
sL=rmsL*s(ind+in);
sR=rmsR*s(ind);
    
sL=ramp(sL,t,rampInd);
sR=ramp(sR,t,rampInd);