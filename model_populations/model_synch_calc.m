group = 150;
fol = 'f:\desktop\WorkingFolder\model_brian';
load([fol '\ICxSpikes_ITD' num2str(group) '.mat'])
spikesICx = cellfun(@(x) x/1000, spikesICx, 'UniformOutput', 0);
% %
Window = sort([-0.1 0.1]);
% binWidth = 0.005; int = 5;
binWidth = 0.001; int = 15;
binNum = length(Window(1):binWidth:Window(2))-2;
pairs = combnk(1:size(spikesICx, 1), 2);
dur = 0.3;

D = cell(length(pairs), size(spikesICx, 3));

for param = 1:size(spikesICx, 3)
    disp(['params = ' num2str(param)])
    for p = 1:length(pairs)
        disp(num2str(p))
        st_u1 = spikesICx(pairs(p, 1), :, param);
        st_u2 = spikesICx(pairs(p, 2), :, param);
        
        sc_u1 = cellfun(@length, st_u1);
        sc_u2 = cellfun(@length, st_u2);

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
        
        holder = cell(1,6);
        holder{1} = [pairs(p,1),pairs(p,2),param];
        for i = 2:3
            holder{i} = nan(1,depvar);
        end
        for i = 4:6
            holder{i}=cell(1,depvar);
        end
        
        for t = 1:depvar
            
            rate1=mean(sc_u1(t,:)/ dur);
            rate2=mean(sc_u2(t,:)/ dur);
            gm = sqrt(rate1*rate2);
            
            SmCCG=(smooth(CCG(t,:))/gm)/binNum;
            SmShiftCCG=(smooth(shiftCCG(t,:))/gm)/binNum;
            
            SmCCGCorrected=SmCCG-SmShiftCCG;
            
            holder{2}(1,t)=rate1;
            holder{3}(1,t)=rate2;
            holder{4}{1,t}=SmCCG;
            holder{5}{1,t}=SmShiftCCG;
            holder{6}{1,t}=SmCCGCorrected;
            
            clear positive pos less more Low High bandPEAK sens1 sens2
            
        end     % for "i" the second time
        
        D{p, param} = holder;
        
        clear holder depvar rate* i PEAK check SmCCG curves CORR trial flank* sumband SmShiftCCG shiftCCG SmCCGCorrected CCG
        
    end
end

synch = D;

save(['synch_model_ITD' num2str(group)], 'synch', 'pairs')