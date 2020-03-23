function power=getPSD(bandwidth, PSD,w)

Nf=length(w); % number of frequencies in signal
power=zeros(1,Nf);

if PSD==1
	power=exp(-w.^2/(2*bandwidth^2));		% Gaussian power spectrum   
elseif PSD==2
   ind=find(abs(w)<=bandwidth);
   power(ind)=1;   
elseif PSD==3
   tRC=bandwidth/(2*pi);
   power=1./(1+w.^2*tRC);
%elseif PSD==4
%   ind=find(abs(w)<=bandwidth & abs(w)>bmin);
%   power(ind)=1;
end

   
