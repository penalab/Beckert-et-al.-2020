% cd('f:\desktop\Data\Analysis\Data')
% cd('Z:\Michael Beckert\data\Analysis\Data')
% load('FreeFieldRawDATA_notsidecorrected_FFffomitted.mat','OT')
%%

cd('f:\desktop\WorkingFolder\distance\data')
load('FreeField_OT_wDelay','OT')
% load('FreeField_OT_wDelay_unsorted','OT')
data = OT; clear OT
cd('F:\Desktop\WorkingFolder')

data = data_trim(data, 0.1);

distanceES = cell(200,1);
distanceEuc = cell(200,1);
distanceVP = cell(200,1);
Mean = cell(200,1);
anID = cell(200,1);

q = 100;
dur = 0.15;
del = 0.05;

start = 0.1;

count = 1;

for date = 1:length(data)
    for site = 1:length(data{date})
        for unit = 2:length(data{date}{site})
            
            d = data{date}{site}{unit};
            disp(['Date - ' num2str(date) '; Site - ' num2str(site) '; Unit - ' num2str(unit)]);
            filetypes = fieldnames(d);
            idx = regexpi(filetypes,'FFf.');
            filetypes = filetypes(~cellfun(@isempty,idx));
            
            distanceEuc{count}=cell([length(filetypes) 1]);
            distanceVP{count}=cell([length(filetypes) 1]);
            distanceES{count}=cell([length(filetypes) 1]);
            Mean{count}=cell([length(filetypes) 1]);
            anID{count} = data{date}{site}{1}(1:3);
            check = 0;
            
            for file=1:length(filetypes)
                TrialData = d.(filetypes{file}).spikes_times;
                
                for i = 1:length(TrialData(:))
                    TrialData{i}(TrialData{i} < start) = [];
                end
                
                distanceEuc{count}{file} = nan(1, size(TrialData,1));
                distanceVP{count}{file} = nan(1, size(TrialData,1));
                distanceES{count}{file} = nan(1, size(TrialData,1));
                Mean{count}{file} = nanmean(cellfun(@length, TrialData), 2)';
                
                for u = 1:size(TrialData, 1)
                    holder = nan(size(TrialData, 2));
                    holder2 = nan(size(TrialData, 2));
                    holder3 = nan(size(TrialData, 2));
                    for t1 = 1:size(TrialData, 2)
                        for t2 = 1:size(TrialData, 2)
                            if ~isempty(TrialData{u, t1}) && ~isempty(TrialData{u, t2})
%                             holder(t1, t2) = euclidean_distance(TrialData{u, t1}, TrialData{u, t2}, start:0.001:del + dur);
                            holder2(t1, t2) = VPdistance_reich(TrialData{u, t1}, TrialData{u, t2}, q);
%                             [holder3(t1, t2), ~] = Event_Sync(TrialData{u, t1}, TrialData{u, t2});
                            end
                        end
                    end
                    idx = ones(size(TrialData, 2));
                    idx = tril(idx, -1);
                    idx = find(idx == 1);
                    distanceEuc{count}{file}(u) = nanmean(holder(idx));
                    distanceVP{count}{file}(u) = nanmean(holder2(idx));
                    distanceES{count}{file}(u) = nanmean(holder3(idx));
                end
                
                clear TrialData trialdata u t1 t2 holder holder2 holder3 idx Indexer
                check = check + sum(sum(Mean{count}{file}));
            end
            
            if check > 0
                count = count + 1;
            end
            
        end
    end
end

L = cellfun(@isempty,distanceES);
distanceEuc(L)=[];
distanceVP(L)=[];
distanceES(L)=[];
Mean(L)=[];
anID(L) = [];

for i = 1:length(distanceES)
    distanceEuc{i} = cell2mat(distanceEuc{i});
    distanceVP{i} = cell2mat(distanceVP{i});
    distanceES{i} = cell2mat(distanceES{i});
    Mean{i} = cell2mat(Mean{i});
end

distanceEuc = cellfun(@(x) nanmean(x,1),distanceEuc,'UniformOutput',0);
distanceVP = cellfun(@(x) nanmean(x,1),distanceVP,'UniformOutput',0);
distanceES = cellfun(@(x) nanmean(x,1),distanceES,'UniformOutput',0);
Mean = cellfun(@(x) nanmean(x,1),Mean,'UniformOutput',0);

clearvars -except Index distanceES distanceEuc distanceVP Mean

% cd('f:\desktop\Data\Analysis\ReproducibilityPaper\working')
cd('f:\desktop\WorkingFolder')

save otFF_distance_2000_omit
% %
cd('f:\desktop\WorkingFolder')
load('f:\desktop\WorkingFolder\tetrodekeys\SpeakerIndex.mat','SpeakerIndex')

for t = 1:2
    load('otFF_distance_2000_omit.mat');
    
    for i = 1:length(distanceES)
        
        best = min(SpeakerIndex(Mean{i} == max(Mean{i}),t));
        idx = SpeakerIndex(:,t) == best;
        
        distanceEuc{i} = distanceEuc{i}(idx)';
        distanceVP{i} = distanceVP{i}(idx)';
        distanceES{i} = distanceES{i}(idx)';
        Mean{i} = Mean{i}(idx)';
        
        clear idx best
        
    end
    
    clear i
    
    switch t
        case 1
            save otAZ_distance_2000_omit distanceES distanceEuc distanceVP Mean
        case 2
            save otEL_distance_2000_omit distanceES distanceEuc distanceVP Mean
    end
    
end

clear t
