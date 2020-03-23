q = 200;
del = 0.05;

for type = 1
    switch type
        case 1
%             cd('f:\desktop\Data\ICls\Data\RawData\FF')
            cd('Z:\Michael Beckert\data\ICls\Data\RawData\FF')
            dur = 0.5;
            start = 0.1;
        case 2
%             cd('f:\desktop\Data\ICls\Data\RawData\Dichotic\ILD')
            cd('Z:\Michael Beckert\data\ICls\Data\RawData\Dichotic\ILD')
            dur = 0.3;
            start = 0.1;
        case 3
%             cd('f:\desktop\Data\ICls\Data\RawData\Dichotic\ITD')
            cd('Z:\Michael Beckert\data\ICls\Data\RawData\Dichotic\ITD')
            dur = 0.3;
            start = 0.1;
    end
    
    FILES = ls;
    FILES = FILES(3:end,:);
    
    distanceEuc = cell(size(FILES, 1), 1);
    distanceVP = cell(size(FILES, 1), 1);
    distanceES = cell(size(FILES, 1), 1);
    Index = cell(size(FILES, 1), 1);
    Mean = cell(size(FILES, 1), 1);
    
    for file=1:size(FILES,1)
        
        if type >1
            load(FILES(file,:), 'curvedata')
            TrialData = curvedata.spike_times;
            clear curvedata
        else
            load(FILES(file,:), 'TrialData', 'Indexer')
            Index{file}=Indexer;
        end

        TrialData = cellfun(@(x) x/1000, TrialData, 'UniformOutput', 0);
        TrialData = cellfun(@transpose, TrialData, 'UniformOutput', 0);
        
        for i = 1:length(TrialData(:))
            TrialData{i}(TrialData{i} < start) = [];
        end
        
        Mean{file} = nanmean(cellfun(@length, TrialData), 2);   
        
        distanceEuc{file} = nan(size(TrialData, 1), 1);
        distanceVP{file} = nan(size(TrialData, 1), 1);
        distanceES{file} = nan(size(TrialData, 1), 1);
        for u = 1:size(TrialData, 1)
            disp(['type = ' num2str(type) ', file = ' num2str(file) ', unit = ' num2str(u)])
            tic
            holder = nan(size(TrialData, 2));
            holder2 = nan(size(TrialData, 2));
            holder3 = nan(size(TrialData, 2));
            
            idx = ones(size(TrialData, 2));
            idx = tril(idx, -1);
            
            for t1 = 1:size(TrialData, 2)
                for t2 = 1:size(TrialData, 2)
                    if idx(t1, t2) == 1
                        if ~isempty(TrialData{u, t1}) && ~isempty(TrialData{u, t2})
                        % holder(t1, t2) = euclidean_distance(TrialData{u, t1}, TrialData{u, t2}, start:0.001:del + dur);
                        holder2(t1, t2) = VPdistance_reich(TrialData{u, t1}, TrialData{u, t2}, q);
                        % [holder3(t1, t2), ~] = Event_Sync(TrialData{u, t1}, TrialData{u, t2});
                        end
                    end
                end
            end
            idx = find(idx == 1);
            distanceEuc{file}(u) = nanmean(holder(idx));
            distanceVP{file}(u) = nanmean(holder2(idx));
            distanceES{file}(u) = nanmean(holder3(idx));
            toc
        end
        
        clear TrialData trialdata u t1 t2 holder idx Indexer holder2 holder3
        
    end
    
    cd('f:\desktop\WorkingFolder')
    
    switch type
        case 1
            save(['iclsFF_distance_' num2str(q) '_omit'], 'distanceEuc', 'distanceVP', 'distanceES', 'Index', 'Mean')
            for t = 1:2
            load(['iclsFF_distance_' num2str(q) '_omit.mat'])
                for i = 1:length(Index)
                    idx = find(Index{i}(t,:) == mode(Index{i}(t,:)));
                    distanceEuc{i} = distanceEuc{i}(idx);
                    distanceVP{i} = distanceVP{i}(idx);
                    distanceES{i} = distanceES{i}(idx);
                    Index{i} = Index{i}(:,idx);
                    Mean{i} = Mean{i}(idx);
                end
                clear i idx
                switch t
                    case 1 
                        save(['iclsAZ_distance_' num2str(q) '_omit'], 'distanceEuc', 'distanceVP', 'distanceES', 'Index', 'Mean')
                    case 2
                        save(['iclsEL_distance_' num2str(q) '_omit'], 'distanceEuc', 'distanceVP', 'distanceES', 'Index', 'Mean')
                end
            end
            clear t
        case 2
            save(['iclsILD_distance_' num2str(q) '_omit'], 'distanceEuc', 'distanceVP', 'distanceES', 'Index', 'Mean')
        case 3
            save(['iclsITD_distance_' num2str(q) '_omit'], 'distanceEuc', 'distanceVP', 'distanceES', 'Index', 'Mean')
    end
    
end

clear type bins FILES file