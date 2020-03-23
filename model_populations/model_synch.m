%%
load('f:\desktop\Code\python\reproducibility\data\ICxSpikes_ITD150_2_highSNR.mat')

Window = sort([-0.1 0.1]);
binWidth = 0.001; int = 15;
binNum = length(Window(1):binWidth:Window(2))-2;

dur = 0.3;

st_u1 = spikesICx1;
st_u2 = spikesICx2;
sc_u1 = cellfun(@length, st_u1);
sc_u2 = cellfun(@length, st_u2);

[depvar, reps] = size(st_u1);

CCG = zeros(depvar,binNum);
shiftCCG = zeros(depvar,binNum);
                
for t = 1:depvar
    for r = 1:reps
        if ~isempty(st_u1{t,r})
            for spike = 1:length(st_u1{t,r})
                times = st_u1{t,r}(spike) - st_u2{t,r};
                his=hist(times,Window(1):binWidth:Window(2)); %% Here is where you can change binwidth
                CCG(t,:)=CCG(t,:)+his(2:end-1);
                
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
                
holder = cell(1,7);
holder{1} = [1, 2];
for i = 2:4
    holder{i} = nan(1,depvar);
end
for i = 5:7
    holder{i}=cell(1,depvar);
end

for t = 1:depvar
    
    rate1=mean(sc_u1(t,:)/dur);
    rate2=mean(sc_u2(t,:)/dur);
    
    SmCCG=(smooth(CCG(t,:))/sqrt(rate1*rate2))/binNum;
    SmShiftCCG=(smooth(shiftCCG(t,:))/sqrt(rate1*rate2))/binNum;
    
    SmCCGCorrected=SmCCG-SmShiftCCG;
    
    flankS=std([SmCCGCorrected(1:round((binNum/5)));SmCCGCorrected(end-round((binNum/5)):end)]);
    flankM=mean([SmCCGCorrected(1:round((binNum/5)));SmCCGCorrected(end-round((binNum/5)):end)]);
    center = SmCCGCorrected(round(binNum/2) - int:round(binNum/2) + int);
    
    PEAK=SmCCGCorrected(round(binNum/2));
    check1 = PEAK>flankM+5*flankS;
    
    holder{2}(1,t)=rate1;
    holder{3}(1,t)=rate2;
    holder{4}(1,t)=check1;
    holder{5}{1,t}=SmCCG;
    holder{6}{1,t}=SmShiftCCG;
    holder{7}{1,t}=SmCCGCorrected;
    
    
end     % for "i" the second time

icx_synch = holder;

