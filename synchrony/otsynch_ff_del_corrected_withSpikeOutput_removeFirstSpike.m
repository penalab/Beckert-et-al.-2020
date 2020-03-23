%%
% load('f:\desktop\Data\Analysis\Data\FreeFieldRawDATA_notsidecorrected_FFffomitted.mat','OT');
load('K:\WorkingFolder\distance\data\FreeField_OT_wDelay','OT')
data = OT; clear OT

load('K:\WorkingFolder\tetrodekeys\blacklist.mat')

Window = sort([-0.1 0.1]);
% binWidth = 0.005; int = 5;
binWidth = 0.001; int = 15;
binNum = length(Window(1):binWidth:Window(2))-2;
norminc = 0:0.1:1;

dur = 0.2;
start = 0.1;
pattern = 'frozen.|b.|i.';                   % These are to remove
pattern2 = '.';                   % These are to keep
f = 'ff';

D=[];

for date=1:size(data,1)
    for site = 1:length(data{date})
        list=combnk(2:length(data{date}{site}),2);
        for pair = 1:size(list,1)
            
            disp(['date: ' num2str(date) '; site: ' num2str(site) '; unit1: ' num2str(list(pair,1)) '; unit2: ' num2str(list(pair,2))]);
            
            unit1 = list(pair,1); unit2 = list(pair,2);
            chck1 = blacklist(:,1) == date;
            chck2 = blacklist(:,2) == site;
            chck3 = blacklist(:,3) == unit1;
            chck4 = blacklist(:,3) == unit2;
            check = chck1 .* chck2 .* (chck3 + chck4);
            
            if sum(check) > 0
                pattern2 = 'FFf1';
                disp('caught')
            elseif sum(check) == 0
                pattern2 = '.';
            end
            
            unit1 = data{date}{site}{list(pair,1)};
            unit2 = data{date}{site}{list(pair,2)};
            filetypes1 = fieldnames(unit1);
            filetypes2 = fieldnames(unit2);
            filetypes = intersect(filetypes1,filetypes2);
            idx = regexpi(filetypes,pattern);
            idx = cellfun(@isempty,idx);
            filetypes = filetypes(idx);
            idx = regexpi(filetypes,pattern2);
            idx = ~cellfun(@isempty,idx);
            filetypes = filetypes(idx);
            
            clear filetypes1 filetypes2 idx chck* check
            
            if ~isempty(filetypes)
                st_u1 = []; st_u2 = []; b_u1 = []; b_u2 = [];
                for file = 1:length(filetypes)
                    st_u1 = [st_u1,unit1.(filetypes{file}).spikes_times];
                    st_u2 = [st_u2,unit2.(filetypes{file}).spikes_times];
                    b_u1 = [b_u1, unit1.(filetypes{file}).base_times];
                    b_u2 = [b_u2, unit1.(filetypes{file}).base_times];
                end
                st_u1_orig = st_u1;
                st_u2_orig = st_u2;
                
                for i = 1:length(st_u1(:))
                    st_u1{i}(st_u1{i} < start) = [];
                    st_u2{i}(st_u2{i} < start) = [];
                end
                
%                 for i = 1:length(st_u1(:))
%                     if ~isempty(st_u1{i})
%                         st_u1{i}(1) = [];
%                     end
%                     if ~isempty(st_u2{i})
%                         st_u2{i}(1) = [];
%                     end
%                 end     
                % Would you also like to remove a second spike? Just to be
                % safe?
%                 for i = 1:length(st_u1(:))
%                     if ~isempty(st_u1{i})
%                         st_u1{i}(1) = [];
%                     end
%                     if ~isempty(st_u2{i})
%                         st_u2{i}(1) = [];
%                     end
%                 end                          
                sc_u1 = cellfun(@length, st_u1);
                sc_u2 = cellfun(@length, st_u2);
                b_u1 = cellfun(@length, b_u1);
                b_u2 = cellfun(@length, b_u2);
                thresh = sqrt(b_u1 .* b_u2);
                thresh = nanmean(thresh(:)) + 2 * nanstd(thresh(:));
                
                [depvar, reps] = size(sc_u1);
                
                CCG = zeros(depvar,binNum);
                shiftCCG = zeros(depvar,binNum);
                
                for t = 1:depvar
                    for r = 1:reps
                        if ~isempty(st_u1{t,r})
                            for spike = 1:length(st_u1{t,r})
                                times = st_u1{t,r}(spike) - st_u2{t,r};
                                his=hist(times,Window(1):binWidth:Window(2)); %% Here is where you can change binwidth
                                CCG(t,:)=CCG(t,:)+his(2:end-1);
                                
                                % regular shift
                                if r==size(sc_u1,2)
                                    times = st_u1{t,r}(spike) - st_u2{1};
                                else
                                    times = st_u1{t,r}(spike) - st_u2{t,r+1};
                                end
                                
                                his=hist(times,Window(1):binWidth:Window(2));
                                shiftCCG(t,:)=shiftCCG(t,:)+his(2:end-1);
                            end
                        end
                    end
                end % for "t"
                
                clear his spike times
                
                holder = cell(1,15);
                holder{1} = [date,site,list(pair,1),list(pair,2)];
                for i = 2:10
                    holder{i} = nan(1,depvar);
                end
                for i = 11:13
                    holder{i}=cell(1,depvar);
                end
                
                for t = 1:depvar
                    
                    rate1=mean(sc_u1(t,:)/ (dur - start));
                    rate2=mean(sc_u2(t,:)/ (dur - start));
                    gm = sqrt(rate1*rate2);
                    
                    SmCCG=(smooth(CCG(t,:))/gm)/binNum;
                    SmShiftCCG=(smooth(shiftCCG(t,:))/gm)/binNum;
                    
                    SmCCGCorrected=SmCCG-SmShiftCCG;
                    
                    flankS=std([SmCCGCorrected(1:round((binNum/5)));SmCCGCorrected(end-round((binNum/5)):end)]);
                    flankM=mean([SmCCGCorrected(1:round((binNum/5)));SmCCGCorrected(end-round((binNum/5)):end)]);
                    center = SmCCGCorrected(round(binNum/2) - int:round(binNum/2) + int);
                    
                    PEAK=SmCCGCorrected(round(binNum/2));
                    check1 = PEAK>flankM+5*flankS;
                    positive=(SmCCGCorrected>=(PEAK/2));
                    pos=diff(positive);
                    less=pos(1:round(binNum/2));
                    more=pos(round(binNum/2):end);
                    Low=find(less==1,1,'last');
                    High=find(more==-1,1,'first');
                    Low=length(less)-Low;
                    bandPEAK=High+Low;
                    
                    if isempty(bandPEAK)
                        bandPEAK=nan;
                        sumband=nan;
                    else
                        sumband=sum(SmCCGCorrected(length(less)-Low:length(less)+High))/bandPEAK;
                    end
                    
                    holder{2}(1,t)=PEAK;
                    holder{3}(1,t)=rate1;
                    holder{4}(1,t)=rate2;
                    holder{5}(1,t)=check1;
                    holder{6}(1,t)=gm > thresh;
                    %                     holder{7}(1,i)=check;
                    holder{8}(1,t)=bandPEAK;
                    holder{9}(1,t)=sum(center);
                    holder{10}(1,t)=sumband;
                    holder{11}{1,t}=SmCCG;
                    holder{12}{1,t}=SmShiftCCG;
                    holder{13}{1,t}=SmCCGCorrected;
                    
                    clear positive pos less more Low High bandPEAK sens1 sens2
                    
                end     % for "i" the second time
                
                holder{14}{1} = st_u1_orig;
                holder{15}{1} = st_u2_orig;
                
                if isfield(D,f)
                    for i = 1:10
                        if size(D.(f){i},2) > size(holder{i},2)
                            padding = nan(1,size(D.(f){i},2)-size(holder{i},2));
                            D.(f){i}=[D.(f){i};holder{i},padding];
                        elseif size(D.(f){i},2) < size(holder{i},2)
                            padding = nan(size(D.(f){i},1),size(holder{i},2)-size(D.(f){i},2));
                            D.(f){i}=[D.(f){i},padding;holder{i}];
                        else
                            D.(f){i}=[D.(f){i};holder{i}];
                        end
                    end
                    for i = 11:13
                        if size(D.(f){i},2) > size(holder{i},2)
                            padding = cell(1,size(D.(f){i},2)-size(holder{i},2));
                            D.(f){i}=[D.(f){i};holder{i},padding];
                        elseif size(D.(f){i},2) < size(holder{i},2)
                            padding = cell(size(D.(f){i},1),size(holder{i},2)-size(D.(f){i},2));
                            D.(f){i}=[D.(f){i},padding;holder{i}];
                        else
                            D.(f){i}=[D.(f){i};holder{i}];
                        end
                    end
                    for i = 14:15
                       D.(f){i} = [D.(f){i};holder{i}]; 
                    end
                else
                    for i = 1:length(holder)
                        D.(f){i}=holder{i};
                    end
                end
                
                clear holder depvar rate* i PEAK check SmCCG curves CORR trial flank* sumband SmShiftCCG shiftCCG SmCCGCorrected CCG
                
            end         % for "file"
            
            clear filetypes
            
        end             % for "pair"
    end                 % for "site"
end                     % for "date"

ot_synch = D;

clearvars -except ot_synch

cd('K:\WorkingFolder\synchony')
save ot_synch_ff_del_fr_200_cor_spiketimes