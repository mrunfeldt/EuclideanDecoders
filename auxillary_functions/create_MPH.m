
% % % Creates Modulation Period Histogram for a single neuron/channel % % % 
% INPUTS:(1) "raw" : spikeTimes in seconds {fM}{trial}
% (2) "fMs" : vector of modulation frequencies (in Hz)
% (3) "duration" : duration of stimulus in seconds
% (4) varargin: set to 1 if you want plots produced, else leave empty.

% %  Updated 2015_04_05: binsize changes as a function of modulation
% frequency (hardcoded). MJRunfeldt 

function [MPH] = create_MPH(raw,fMs,duration,varargin)

% input = raw ; % {fM}{tr}
% fMs = data.fM ; 

nBinPerCycle = 80 ; % temporal resolution of MPH (relates to bin sizes)
binSizes = (1./fMs) / nBinPerCycle ; % customize binSize per 

MPH = cell(1,length(fMs));
for f = 1:length(fMs)
    nPerRec = fMs(f) / duration ; % number of cycles in recording
    period = 1 / fMs(f) ; % duration of period in seconds
    allSpikes = cell2mat(raw{f}) ; % all spike times all trials
    nTrials = length(raw{f});
    
    if  nPerRec < 1e3 % %  limit so it doesn't run forever
    on = 0;nSpikes=[];
    for cy = 1:nPerRec
        chunk= on:binSizes(f):on+period ;on=on+period;
        dum = histc(allSpikes,chunk) ; dum(end)=[];
        nSpikes(cy,:) = dum;
    end
    MPH{f} = sum(nSpikes) ./ cy;
    
    else MPH{f} = 0;
    end
    %figure;plot(MPH{f});title(num2str(fMs(f)))
end



if nargin >3
figure;nsubsX = ceil(f/4); nsubsY = ceil(f/nsubsX) ;
    for ff = 1:f
subplot(nsubsX,nsubsY,ff);
plot(MPH{ff});title([num2str(fMs(ff)), ' Hz'] )
p1=length(MPH{ff}); xlim([0 p1]);per = round(1e4/fMs(ff))/10 ;
set(gca,'xtick',[0 p1],'xticklabel',{'0',num2str(per)})
    end

end