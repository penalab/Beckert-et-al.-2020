cd('K:\WorkingFolder')

% % run correction for icls
% 1 - ITD
% 2 - ILD
% 3 - AZ
% 4 - EL

icls = cell(4, 1);

for t = 1:size(icls, 1)
    switch t
        case 1
            load('K:\WorkingFolder\distance\altcorrections\distance_correction_0_p4_c020_q_100min.mat', 'key', 'cfun')
            load('K:\WorkingFolder\distance\varyingQ\iclsITD_distance_100_omit.mat', 'distanceEuc', 'distanceVP', 'distanceES', 'Mean')
        case 2
            load('K:\WorkingFolder\distance\altcorrections\distance_correction_0_p4_c020_q_100min.mat', 'key', 'cfun')
            load('K:\WorkingFolder\distance\varyingQ\iclsILD_distance_100_omit.mat', 'distanceEuc', 'distanceVP', 'distanceES', 'Mean')
        case 3
            load('K:\WorkingFolder\distance\altcorrections\distance_correction_0_p4_c040_q_100min.mat', 'key', 'cfun')
            load('K:\WorkingFolder\distance\varyingQ\iclsAZ_distance_100_omit.mat', 'distanceEuc', 'distanceVP', 'distanceES', 'Mean')
        case 4
            load('K:\WorkingFolder\distance\altcorrections\distance_correction_0_p4_c040_q_100min.mat', 'key', 'cfun')
            load('K:\WorkingFolder\distance\varyingQ\iclsEL_distance_100_omit.mat', 'distanceEuc', 'distanceVP', 'distanceES', 'Mean')
    end
    
    cor_vp = @(x) cfun{1}(1) * x .^ 4 + cfun{1}(2) * x .^ 3 + cfun{1}(3) * x .^ 2 + cfun{1}(4) * x + cfun{1}(5);
    cor_euc = @(x) cfun{2}(1) * x .^ 4 + cfun{2}(2) * x .^ 3 + cfun{2}(3) * x .^ 2 + cfun{2}(4) * x + cfun{1}(5);
    cor_es = @(x) cfun{3}(1) * x .^ 4 + cfun{3}(2) * x .^ 3 + cfun{3}(3) * x .^ 2 + cfun{3}(4) * x + cfun{1}(5);

    es = cell(size(distanceES));
    euc = cell(size(distanceEuc));
    vp = cell(size(distanceVP));
    
    for u = 1:length(Mean)
        es{u} = nan(size(distanceES{u}));
        euc{u} = nan(size(distanceEuc{u}));
        vp{u} = nan(size(distanceVP{u}));
        for i = 1:size(Mean{u}, 1)
            
            es{u}(i) = distanceES{u}(i) - cor_es(Mean{u}(i));
            euc{u}(i) = distanceEuc{u}(i) - cor_euc(Mean{u}(i));
%             if Mean{u}(i) >= 1
                vp{u}(i) = distanceVP{u}(i) - cor_vp(Mean{u}(i));
%             elseif Mean{u}(i) < 1
%                 vp{u}(i) = nan;
%             end
            
        end
    end
    
    icls{t}.euc = euc;
    icls{t}.vp = vp;
    icls{t}.es = es;
    icls{t}.rf = Mean;
    
end

clear distanceEuc distanceVP distanceES Mean t u i idx es euc vp

% % Do correction for ot data
% 1 - ITD (unfrozen)
% 2 - ILD (unfrozen)
% 3 - AZ
% 4 - EL
% 5 - froz (dichotic frozen)
% 6 - array (dichotic frozen)

ot = cell(6, 1);

for t = 1:size(ot, 1)
    switch t
        case 1
            load('K:\WorkingFolder\distance\altcorrections\distance_correction_0_p4_c010_q_100.mat', 'key', 'cfun')
            load('K:\WorkingFolder\distance\data\distance_OTdichotic_itd_noonset.mat', 'euc_c', 'vp_c', 'es_c', 'rm_c')
            E = euc_c; V = vp_c; S = es_c; M = rm_c;
            clear ITD euc_c vp_c es_c rm_c
        case 2
            load('K:\WorkingFolder\distance\altcorrections\distance_correction_0_p4_c010_q_100.mat', 'key', 'cfun')
            load('K:\WorkingFolder\distance\data\distance_OTdichotic_ild_noonset.mat', 'euc_c', 'vp_c', 'es_c', 'rm_c')
            E = euc_c; V = vp_c; S = es_c; M = rm_c;
            clear ILD euc_c vp_c es_c rm_c
        case 3
            load('K:\WorkingFolder\distance\altcorrections\distance_correction_0_p4_c010_q_100min.mat', 'key', 'cfun')
            load('K:\WorkingFolder\distance\varyingQ\otAZ_distance_100_omit.mat', 'distanceEuc', 'distanceVP', 'distanceES', 'Mean')
            E = distanceEuc; V = distanceVP; S = distanceES; M = Mean;
            clear distanceEuc distanceVP distanceES Mean
        case 4
            load('K:\WorkingFolder\distance\altcorrections\distance_correction_0_p4_c010_q_100min.mat', 'key', 'cfun')
            load('K:\WorkingFolder\distance\varyingQ\otEL_distance_100_omit.mat', 'distanceEuc', 'distanceVP', 'distanceES', 'Mean')
            E = distanceEuc; V = distanceVP; S = distanceES; M = Mean;
            clear distanceEuc distanceVP distanceES Mean
        case 5
            load('K:\WorkingFolder\distance\altcorrections\distance_correction_0_p4_c020_q_100.mat', 'key', 'cfun')
            load('K:\WorkingFolder\distance\data\distance_OTdichotic_froz_noonset.mat', 'euc_c', 'vp_c', 'es_c', 'rm_c')
            E = euc_c; V = vp_c; S = es_c; M = rm_c;
            clear FROZ euc_c vp_c es_c rm_c
        case 6
            load('K:\WorkingFolder\distance\altcorrections\distance_correction_0_p4_c020_q_100.mat', 'key', 'cfun')
            array = load('K:\WorkingFolder\distance\data\arrayITD_distance_noonset.mat');
            idx = cellfun(@(x) length(size(x)) ~= 2, array.rf);
            array.es(idx) = []; array.euc(idx) = []; array.vp(idx) = []; array.rf(idx) = [];
            idx = cellfun(@isempty, array.rf);
            array.es(idx) = []; array.euc(idx) = []; array.vp(idx) = []; array.rf(idx) = [];
            array.rf = cellfun(@transpose, array.rf, 'UniformOutput', 0);
            E = array.euc; V = array.vp; S = array.es; M = array.rf;
            clear array idx
    end
    
    cor_vp = @(x) cfun{1}(1) * x .^ 4 + cfun{1}(2) * x .^ 3 + cfun{1}(3) * x .^ 2 + cfun{1}(4) * x + cfun{1}(5);
    cor_euc = @(x) cfun{2}(1) * x .^ 4 + cfun{2}(2) * x .^ 3 + cfun{2}(3) * x .^ 2 + cfun{2}(4) * x + cfun{1}(5);
    cor_es = @(x) cfun{3}(1) * x .^ 4 + cfun{3}(2) * x .^ 3 + cfun{3}(3) * x .^ 2 + cfun{3}(4) * x + cfun{1}(5);

    % Apply the correction
    
    for u = 1:length(M)
        for s = 1:size(M{u}, 1)
            
            S{u}(s) = S{u}(s) - cor_es(M{u}(s));
            E{u}(s) = E{u}(s) - cor_euc(M{u}(s));
%             if M{u}(s) >= 1
                V{u}(s) = V{u}(s) - cor_vp(M{u}(s));
%             elseif M{u}(s) < 1
%                 V{u}(s) = nan;
%             end
            
        end
    end
    
    ot{t}.euc = E;
    ot{t}.vp = V;
    ot{t}.es = S;
    ot{t}.rf = M;
    
end

clear t u s E V S M

% % Run correction for the special case data
% 1 - for individual cells across speakers
% 2 - for recording sites across neurons

% special = cell(2, 1);
% 
% load('f:\desktop\WorkingFolder\distance\data\distance_correction_c_01_10_long.mat', 'key', 'cfun')
% 
% cor_vp = @(x) cfun{1}(1) * x + cfun{1}(2);
% cor_euc = @(x) cfun{2}(1) * x + cfun{2}(2);
% cor_es = @(x) cfun{3}(1) * x + cfun{3}(2);
% 
% for t = 1:size(special, 1)
%     switch t
%         case 1
%             load('f:\desktop\WorkingFolder\distance\data\distance_ot_betweenspeakers_noonset.mat', 'distanceEuc', 'distanceVP', 'distanceES', 'Mean')
%             E = distanceEuc; V = distanceVP; S = distanceES; M = Mean;
%             E = cellfun(@transpose, E, 'UniformOutput', 0);
%             V = cellfun(@transpose, V, 'UniformOutput', 0);
%             S = cellfun(@transpose, S, 'UniformOutput', 0);
%             M = cellfun(@transpose, M, 'UniformOutput', 0);
%             clear distanceEuc distanceVP distanceES Mean
%         case 2
%             load('f:\desktop\WorkingFolder\distance\data\distance_betweenneurons_del_noonset.mat', 'EUC', 'VP', 'ES', 'FR1', 'FR2')
%             
%             E = cell(length(VP.FFf1), 1);
%             fr1 = cell(length(VP.FFf1), 1);
%             fr2 = cell(length(VP.FFf1), 1);
%             V = cell(length(VP.FFf1), 1);
%             S = cell(length(VP.FFf1), 1);
%             
%             for p = 1:length(VP.FFf1)
%                 E{p} = [EUC.FFf1{p} , EUC.FFf2{p}];
%                 S{p} = [ES.FFf1{p} , ES.FFf2{p}];
%                 fr1{p} = [FR1.FFf1{p} , FR1.FFf2{p}];
%                 fr2{p} = [FR2.FFf1{p} , FR2.FFf2{p}];
%                 V{p} = [VP.FFf1{p} , VP.FFf2{p}];
%             end
%             
%             E = cellfun(@(x) mean(x, 2), E, 'UniformOutput', 0);
%             S = cellfun(@(x) mean(x, 2), S, 'UniformOutput', 0);
%             fr1 = cellfun(@(x) mean(x, 2), fr1, 'UniformOutput', 0);
%             fr2 = cellfun(@(x) mean(x, 2), fr2, 'UniformOutput', 0);
%             V = cellfun(@(x) mean(x, 2), V, 'UniformOutput', 0);
%                         gm = cellfun(@(x1, x2) sqrt( x1 .* x2 ), fr1, fr2, 'UniformOutput', 0);
%             M = cellfun(@(x1, x2) ( x1 + x2 ) / 2, fr1, fr2, 'UniformOutput', 0);
%             
%             idx = cellfun(@isempty, M);
%             M(idx) = []; 
%             E(idx) = []; S(idx) = []; 
%             V(idx) = [];
%             
%             clear fr1 fr2 EUC ES FR1 FR2 VP p idx
%     end
%     
%     for u = 1:length(M)
%         for s = 1:size(M{u}, 1)
%             
%             S{u}(s) = S{u}(s) - cor_es(M{u}(s));
%             E{u}(s) = E{u}(s) - cor_euc(M{u}(s));
% %             if M{u}(s) >= 1
%                 V{u}(s) = V{u}(s) - cor_vp(M{u}(s));
% %             elseif M{u}(s) < 1
% %                 V{u}(s) = nan;
% %             end
%             
%         end
%     end
%     
%     special{t}.euc = E;
%     special{t}.vp = V;
%     special{t}.es = S;
%     special{t}.rf = M;
%     
% end
% 
% clear t u s E V S M

% %
clear cor* key cfun

% save distance_notcorrected_data2