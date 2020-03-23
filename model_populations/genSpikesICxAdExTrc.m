function [Vm,spikecnt,sp,thresh] = genSpikesICxAdExTrc(ge,gi,I,dt,NeuronParms)

%For one neuron

%% Number of neurons

Lt=length(ge);

%% Time

tt=0:dt:(Lt-1)*dt;

%% Neuron parameters

Tref = NeuronParms(1);
Vrest = NeuronParms(2);
Vreset = NeuronParms(3);
EL = NeuronParms(4);
tRC = NeuronParms(5);
%C = NeuronParms(6);
Ee = NeuronParms(7);
Ei = NeuronParms(8);
th_range = NeuronParms(9);
Del = NeuronParms(10);


tau_th = NeuronParms(11);
alpha = NeuronParms(12);
VT = NeuronParms(13);
ka = NeuronParms(14);
ki = NeuronParms(15);
Vi = NeuronParms(16);

MaxFR=1/Tref; %Maximum firing rate

%eps = dt/C;
eps = dt/tRC;



%% Initialize

inref=0;
spikecnt=0;
spiketimes = zeros(1,floor(dt*Lt*MaxFR)+2); % max possible # spikes.

V=Vrest+rand;

%Record potential for plotting
Vm=ones(1,Lt)*Vreset;
Vm(1) = V;

Vthres = alpha*(V - Vi) + VT + ka*log(1 + exp((V - Vi)/ki));

thresh = zeros(1,Lt);

%% integrate

for n = 2:Lt
    
    %Update threshold   
    th_inf = alpha*(V - Vi) + VT + ka*log(1 + exp((V - Vi)/ki));
    
    Vthres = Vthres + dt/tau_th*(th_inf - Vthres);
    
    %If out of refractory period, update membrane potential
    if inref <= 0
        
        dV=((EL-V)+Del*exp((V-VT)/Del)+I(n-1)+ge(n-1)*(Ee-V)+gi(n-1)*(Ei-V));
        
        V = V + eps*dV;
        
        %Record potential for plotting
        Vm(n)=V;
        
        %Get threshold
        th=Vthres+(2*rand-1)*.5*th_range;
        
        %Test for spikes
        if V > th
            spikecnt = spikecnt + 1;
            
            spiketimes(spikecnt)= tt(n);
            
            inref = Tref; 
            
             %% Spike shape
            %X = [0, Tref/2, Tref];
            %Y = [Vreset, -5, V];
            %p = polyfit(X,Y,2);
            
            V = Vreset;
            %V = th - Dthres;
            
            Vm(n) = -10; %Artificial height of AP 
            
        end
    else
        inref = inref - dt;
        
        %Vm(n) = polyval(p,inref);
    end
    
    thresh(n) = Vthres;
end

%%
sp = spiketimes(1:spikecnt);









