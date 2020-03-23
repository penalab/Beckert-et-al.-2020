
cd('K:/WorkingFolder')
pathstartup_pc

%%

cd('K:\python\model_icx_codeonly\SynchResults')

% % Iterate across variables

% itds = -200:20:-140;
% itds = -160:20:-100;
% itds = -120:20:-60;
% itds = -40:20:20;
% itds = [40:20:80 120];
% itds = 140:20:200;
itds = -200:20:200;
% itds = 160;
repet = 1;

snr_iccl = 2;
snr_icx = 10;
icx_weight = 2;
icx_side = -0.5;
icx_width = -0.5;
freq_width = -0.5;
bw = 0.1;
gainL = 1;
gainU = 1;

fmin = 0.25;
wweight = 11;
bin = 25;

vt = -57;
vi = -67;


for SNR_ICcl = snr_iccl
for SNR_ICx = snr_icx
for ICx_weight = icx_weight
    
for ICx_side = icx_side
for ICx_width = icx_width
for FREQ_width = freq_width
    
for BW = bw
    
for GAINL = gainL
for GAINU = gainU
for W_weight = wweight

        disp('Setting parameters :')
        disp(['SNR_ICcl ' num2str(SNR_ICcl)])
        disp(['SNR_ICx ' num2str(SNR_ICx)])
        disp(['ICx_weight ' num2str(ICx_weight)])
        disp(['ICx_side ' num2str(ICx_side)])
        disp(['ICx_width ' num2str(ICx_width)])
        disp(['freq width ' num2str(FREQ_width)])
        disp(['BW ' num2str(BW)])
        disp(['gainL ' num2str(GAINL)])
        disp(['gainU ' num2str(GAINU)])
        disp(['W_weight ' num2str(W_weight)])
        
        Synch1_Model_Parameters_pipeline(SNR_ICcl, SNR_ICx, ICx_weight, ICx_side, ICx_width, FREQ_width, BW)
        
for rep = repet
for ITD = itds
    
for FMIN = fmin
for BIN = bin
  
for VT = vt
for Vi = vi
    
    disp(['bin ' num2str(BIN)])
    disp(['fmin ' num2str(FMIN)])
    disp(['VT ' num2str(VT)])
    disp(['Vi ' num2str(Vi)])
    
    r_txt = sprintf('%03d', rep);
    
%     if ~isfile(['K:\python\model_icx_codeonly\SynchResults\data\out_synch\modelout_ITD' num2str(ITD) '_' r_txt ...
%     '_SNRICcl' num2str(SNR_ICcl * 10) ...
%     '_SNRICx' num2str(SNR_ICx * 100000) ...
%     '_ICxweight' num2str(ICx_weight * 10) ...
%     '_ICxside' num2str(ICx_side * 10) ...
%     '_ICxwidth' num2str(ICx_width * 10) ...
%     '_freqwidth' num2str(freq_width * 10) ...
%     '_bw' num2str(BW * 100000) ...
%     '_gainL' num2str(GAINL * 100) ...
%     '_gainU' num2str(GAINU) ...
%     '_fmin' num2str(fmin * 100) ...
%     '_Wweight' num2str(W_weight * 10) ...
%     '_bin' num2str(BIN) ...
%     '_VT' num2str(VT) ...
%     '_Vi' num2str(Vi) ...
%     '_shiftshuff_syncalc.mat'])

    disp(['running ' num2str(ITD)])
    params_name = ['K:\python\model_icx_codeonly\SynchResults\data\' ...
        'SynchParms_SNRICcl' num2str(SNR_ICcl * 10) ...
        '_SNRICx' num2str(SNR_ICx * 100000) ...
        '_ICxweight' num2str(ICx_weight * 10) ...
        '_ICxside' num2str(ICx_side * 10) ...
        '_ICxwidth' num2str(ICx_width * 10) ...
        '_freqwidth' num2str(freq_width * 10) ...
        '_bw' num2str(BW * 100000)];
    
%         disp(['calculating cross-correlation : rep = ' num2str(rep) '; ITD = ' num2str(ITD)])
%         Synch2_get_CC_Response(rep, ITD, 'K:\python\model_icx_codeonly\')

%         disp(['calculating ICcl responses : rep = ' num2str(rep) '; ITD = ' num2str(ITD)])
%         Synch3_get_ICcl_Response_reducedSNR(rep, ITD, 'K:\python\model_icx_codeonly\', params_name)
        
        disp('calculating ICx rsesponses')
        disp(['rep = ' num2str(rep) '; ITD = ' num2str(ITD)])
        Synch4_get_ICx_Response_dynamicgain_dynamicSNRICx_pipeline(rep, ITD, 'K:\python\model_icx_codeonly\', SNR_ICcl, SNR_ICx, ICx_weight, ICx_side, ICx_width, FREQ_width, BW, GAINL, GAINU)
        disp('calculating ICx synchrony')
        Synch4b_adjust_ICx_Response_fromshift_RunNCalc_pipeline(rep, ITD, 'K:\python\model_icx_codeonly\', SNR_ICcl, SNR_ICx, ICx_weight, ICx_side, ICx_width, FREQ_width, BW, GAINL, GAINU)
        disp('calculating model output')
        Synch5_Decode_ICxvShiftTimes_trials_pipeline_updated(rep, ITD, 'K:\python\model_icx_codeonly\', SNR_ICcl, SNR_ICx, ICx_weight, ICx_side, ICx_width, FREQ_width, BW, GAINL, GAINU, FMIN, W_weight, BIN, VT, Vi)
        
        var_filename = ['modelout_ITD' num2str(ITD) '_' r_txt ...
            '_SNRICcl' num2str(SNR_ICcl * 10) ...
            '_SNRICx' num2str(SNR_ICx * 100000) ...
            '_ICxweight' num2str(ICx_weight * 10) ...
            '_ICxside' num2str(ICx_side * 10) ...
            '_ICxwidth' num2str(ICx_width * 10) ...
            '_freqwidth' num2str(freq_width * 10) ...
            '_bw' num2str(BW * 100000) ...
            '_gainL' num2str(gainL * 100) ...
            '_gainU' num2str(gainU) ...
            '_fmin' num2str(fmin * 100) ...
            '_Wweight' num2str(W_weight * 10) ...
            '_bin' num2str(bin) ...
            '_VT' num2str(VT) ...
            '_Vi' num2str(Vi) ...
            '_shiftshuff'];

%         disp('calculating model output synchrony')
        Synch5c_calcSynchShift('K:\python\model_icx_codeonly\', var_filename)

        disp('complete')
        disp('----------')


%     end
% disp(['modelout_ITD' num2str(ITD) '_' r_txt ...
%     '_SNRICcl' num2str(SNR_ICcl * 10) ...
%     '_SNRICx' num2str(SNR_ICx * 100000) ...
%     '_ICxweight' num2str(ICx_weight * 10) ...
%     '_ICxside' num2str(ICx_side * 10) ...
%     '_ICxwidth' num2str(ICx_width * 10) ...
%     '_freqwidth' num2str(freq_width * 10) ...
%     '_bw' num2str(BW * 100000) ...
%     '_gainL' num2str(GAINL * 100) ...
%     '_gainU' num2str(GAINU) ...
%     '_shiftshuff.mat'])  

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
end
end

