function [spikes,R,index] = getICclAdExSpikesTrials(X,BITD,Ntri,Ts2,T,SNR,BW,s,NP)

%[spikes,R,index] = getICclAdExSpikesTrials(X,BITD,Ntri,Ts2,T,SNR,BW,s,NP)

%% Parameters

[Nf,Nd,Nt] = size(X);

Nn = length(BITD);

%% Internal ITDs

ITDin=linspace(-200,200,Nd);

%% Find indices of neurons

index = zeros(1,Nn);
for n = 1:Nn
    [~,index(n)] = min(abs(ITDin - BITD(n)));
end
    
%% Initialize rates and spikes

R = zeros(Nf,Nn,Ntri);

spikes = cell(Nf,Nn,Ntri);

%% Get ICcl spikes

for k = 1:Nf

    for n = 1:Nn
     
            
            %Get xcorr rate
            r = zeros(1,Nt);r(:,:) = X(k,index(n),:);
            
            %Transformation
            g = exp(r);
            
            th = mean(g);
        
            g = (g - th);

            g(g<0) = 0;
            
            %RMS of noise signal
            rms_noise = sqrt(rms(g))/SNR;
            
            for i = 1:Ntri
                
                %Generate noise. I had a matching problem, so I generate a
                %longer noise and take the part that matches the length of
                %the signal
                noise = genSignal(T+5,Ts2,rms_noise,2,2*pi*BW);
                noise = noise(1:Nt);

                %Generate the spikes
                [~,fr,sp] = genSpikesICxAdExTrc((g+noise)*s,g*0,g*0,Ts2,NP);

                spikes{k,n,i} = sp(1:fr);

                R(k,n,i) = fr;
                
            end
            
    
    end

end

