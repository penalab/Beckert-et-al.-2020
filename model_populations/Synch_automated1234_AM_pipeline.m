cd('K:\python\model_icx_codeonly\SynchResults')

% % Iterate across variables

% itds = -200:20:200;
itds = 100;
repet = 1;

snr_iccl = 2;
snr_icx = 2;
icx_weight = 1;
icx_side = -0.5;
icx_width = -0.5;
freq_width = -0.5;

for SNR_ICcl = snr_iccl
for SNR_ICx = snr_icx
for ICx_weight = icx_weight
    
for ICx_side = icx_side
for ICx_width = icx_width
for FREQ_width = freq_width
    
        disp('Setting parameters :')
        disp(['SNR_ICcl ' num2str(SNR_ICcl)])
        disp(['SNR_ICx ' num2str(SNR_ICx)])
        disp(['ICx_weight ' num2str(ICx_weight)])
        disp(['ICx_side ' num2str(ICx_side)])
        disp(['ICx_width ' num2str(ICx_width)])
        disp(['freq width ' num2str(FREQ_width)])
        Synch1_Model_Parameters_pipeline(SNR_ICcl, SNR_ICx, ICx_weight, ICx_side, ICx_width, FREQ_width)
        
for rep = repet
for ITD = itds
    
    r_txt = sprintf('%03d', rep);
    
    if ~isfile(['K:\python\model_icx_codeonly\SynchResults\data\testingvars\modelout_ITD' num2str(ITD) '_' r_txt ...
    '_SNRICcl' num2str(SNR_ICcl * 10) ...
    '_SNRICx' num2str(SNR_ICx * 10) ...
    '_ICxweight' num2str(ICx_weight * 10) ...
    '_ICxside' num2str(ICx_side * 10) ...
    '_ICxwidth' num2str(ICx_width * 10) ...
    '_freqwidth' num2str(freq_width * 10) ...
    '_AM_shiftshuff.mat'])
%         disp(['calculating cross-correlation : rep = ' num2str(rep) '; ITD = ' num2str(ITD)])
%         Synch2_get_CC_Response(rep, ITD, 'K:\python\model_icx_codeonly\')

%         disp(['calculating ICcl responses : rep = ' num2str(rep) '; ITD = ' num2str(ITD)])
%         Synch3_get_ICcl_Response_reducedSNR(rep, ITD, 'K:\python\model_icx_codeonly\', SNR_ICcl)
        
        disp('calculating ICx responses')
        disp(['rep = ' num2str(rep) '; ITD = ' num2str(ITD)])
        Synch4b_AM_get_ICx_Response_pipeline(rep, ITD, 'K:\python\model_icx_codeonly\', SNR_ICcl, SNR_ICx, ICx_weight, ICx_side, ICx_width, FREQ_width)
        disp('calculating ICx synchrony')
        Synch4b_AM_adjust_ICx_Response_fromshift_RunNCalc_pipeline(rep, ITD, 'K:\python\model_icx_codeonly\', SNR_ICcl, SNR_ICx, ICx_weight, ICx_side, ICx_width, FREQ_width)
        disp('calculating model output')
        Synch5_Decode_AM_ICxvShiftTimes_trials_pipeline(rep, ITD, 'K:\python\model_icx_codeonly\', SNR_ICcl, SNR_ICx, ICx_weight, ICx_side, ICx_width, FREQ_width)
        disp('complete')
        disp('----------')
    end
end
end

end
end
end

end
end
end

