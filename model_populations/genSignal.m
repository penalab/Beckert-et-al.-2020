function s=genSignal(T,dt,rms,PSD,bandwidth)%% T; total time%% dt; time step%% rms signal level%% Choice of power spectrum type defined by the subroutine getPSD.m%% bandwidth of the signal%% Modified from a version by Brandon Westover%% C. H. Anderson, Sept. 25, 2001 %% %% Copyright (C) by Charles. H. Anderson (All Rights Reserved)%% Dept. Anatomy and Neurobiology%% Washington Univ. School of Medicine%% St. Louis, MO%% cha@shifter.wustl.edu
Nt=2*floor(T/dt/2)+1;			% Number of time samples of s
Nf=Nt;							% # samples doesn't change when converting to freq domain!
df=1/T; 						% Hz
f=[-(Nt-1)/2:(Nt-1)/2]*df; 		% Hz
w=2*pi*f;					    % rad/sec
S=zeros(1,Nt);
power=getPSD(bandwidth, PSD,w);
S=sqrt(power/2).*(randn(1,Nt)+i*randn(1,Nt)); % generate the 'Amplitudes' (signal in the Fourier domain)
s=real(ifft(ifftshift(S)));		% transform to time domain
s=s-mean(s);					% set the mean to sero. (Note that this changes the power spectrum-cha)
s=s/sqrt(mean(s.^2))*rms;		% set the rms