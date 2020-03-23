function [ge_icx,gi_icx] = getICxConductance(spikes,We,Wi,tau_syn,t)

%Get sizes
[Nf,Nn,Nmult] = size(spikes); 
Nt = length(t);
NumN = size(We,1);

%% Get E, I input conductance to ICx neuron

ge_icx = zeros(NumN,Nt);
gi_icx = zeros(NumN,Nt);

for k = 1:Nf
        
        %Initialize
        g_iccl = zeros(Nn,Nt);
        
        for n = 1:Nn
            
                
                %Neurons at same BF, BITD are added
                for i = 1:Nmult
            
                    sp = spikes{k,n,i};
                    Nsp = length(sp);

                    %For each spike, add an exponential 
                    for l = 1:Nsp
                        
                        h = exp(-(t - sp(l))/tau_syn).*(t > sp(l));
                        g_iccl(n,:) = g_iccl(n,:) + h;
                                                
                    end
                end
                        
        end
        
        
        ge_icx = ge_icx + squeeze(We(:,k,:))*g_iccl;
        
        gi_icx = gi_icx + squeeze(Wi(:,k,:))*g_iccl;
            
            
end
