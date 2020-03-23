%% Personal start up command
cd('K:/WorkingFolder')
pathstartup_pc

%% Set list of parameters to put through pipeline

cd('K:\python\model_icx_codeonly\SynchResults')

% % Iterate across variables

% itds = -200:20:-140;

% Parameters for params portion - step 1
% Most are defined 

% Parameters for CC portion - step 2
% None yet

% Parameters for ICcl portion - step 3
snr_iccl = 2;

% Parameters for ICx portion - step 4
snr_icx = 20;
bw = 0.1;
icx_weight = 0.5;
% Keep constant to match general model
    icx_side = -0.5;
    icx_width = -0.5;
    freq_width = -0.5;

% Parameters for Decoding portion - step 5
gainL = 1;
gainU = 1;
fmin = 1;
wweight = 30:5:50;
bin = 20:10:50;


for ITD = itds
for rep = repet
    
    SynchPipeline_step1(
    
for SNR_ICcl = snr_iccl
for SNR_ICx = snr_icx
for ICx_weight = icx_weight
    
for ICx_side = icx_side
for ICx_width = icx_width
for FREQ_width = freq_width
    
for FMIN = fmin
for BW = bw
for GAINL = gainL
for GAINU = gainU
for W_weight = wweight
    
    
    for ITD = itds
    for rep = repet
        
        
    
    
    
    
    end
    end
    
    
end
end
end
end
end

end
end
end


end
end
end

end


