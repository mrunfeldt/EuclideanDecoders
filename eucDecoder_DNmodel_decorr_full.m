
% Bin spike times at high temporal resolution and then convolves binary signal
% with decaying exponential to simulate neural EPSPs.
% Designed for multineuronal convergent model.
% % Applies threshold to EPSPs to generate downstream neuron model; includes (hardcoded) refractory period.
% % Decoding is based on binned raster generated from EPSP model

% % Remove correlations across ALL neurons % % 

% INPUT: (1) spikeTimes {neuron}{stimulus}{trial}
% (2) tau (~full width of EPSP @ half height)
% (3) ISI rebound (4) EPSP threshold (will be converted to fnc of EPSP amp)
% (4) small binSize for EPSP
% (5) binSize for output raster (6) duration of stimulus 
% (7) Number of shuffles for decorr
% % Time unit for all of input must be same (just use seconds)

% MJRunfeldt 2015_08_14: Updated EPSP model to use single eponential,
% Threshold defined relative to maximum amplitude of EPSP, which is:
% Tau / 6.6589
% 2015_09_02: optimized for speed

function [performance, conf_all, spykes,spikeTrains] = ...
    eucDecoder_DNmodel_decorr_full(data,tau,rebound,threshold,...
    binSize_sm,binSize_ras,duration,nShuffles)

% % % % % PARAMETERS % % % % % % % % % % % % % % % % % % % % % % 
isiLim = round(rebound / binSize_sm) ; % ISI limit in small bin samples (HARDCODED)
tau = tau / binSize_sm;
time = [0:round(duration/binSize_sm)]; % time for convolved signal (need integer bin size)
eFun = @(t,tau) t.*(exp(-t.*2.45/tau)); % EPSP function
myExp = eFun(time,tau); % generate exponential
% % % truncate exponential decay signal (otherwise you have a lot of unneeded zeros) % % %
trunc=find(myExp < 1e-4);trunc=trunc(1); if trunc > 10;myExp=myExp(1:trunc);end
rasterTime = [0:binSize_ras:duration] ; % convert to ms for OUTPUT
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % DEFINE THRESHOLD RELATIVE TO MAX AMPLITUDE % % %
maxAmp = tau / 6.6589 ; thresh = threshold*maxAmp ; 

nFMs = length(data{1}) ; % number of modulation frequencies
% spykes = cell(1,nShuffles); spikeTrains = spykes;
% templates = cell(1,nShuffles);
performance = zeros(1,nShuffles);
for f = 1:nFMs % for each mod freq
    nTrials = length(data{1}{f});
    
    for s = 1:nShuffles
% % % Sum NEURONS IN TIME to make convergent mode (use small binSize) % % %    
        merged = zeros(nTrials,length(time)-1); 
    for nn = 1:length(data) % for each neuron
        singleC = data{nn}{f}(randperm(nTrials)) ; % MIX TRIALS
        raz = spikeTimes_toRasters(singleC,binSize_sm,duration) ; % Bin Spike Times SECONDS
        if isempty(raz); raz = zeros(1,length(time)-1);end
        merged = merged+raz; % nTrials x Time
    end % END (nn)
    
% % Generate EPSP model and apply threshold to generate "donwstream neuron" spiketrain % % %    
    boop=[]; % epsp vector (i.e. binary signal convolved with decaying exponential)
        
    for tr=1:nTrials
        epspDum = conv(merged(tr,:),myExp); epspDum = epspDum(1:size(merged,2));
        binary = double(epspDum(1:size(merged,2))>thresh) ;% binary spiketrain
        sTimes = find(binary==1); binDum = binary ; 
        newTimes = sTimes; cnt=1;
% % march forward through binary and set "isiLim" to 0 after each spike
        while cnt < length(newTimes) 
           hit = newTimes(cnt) ; rebound = hit+1 : hit+1+isiLim ; % time frame of rebound
           binDum(rebound)=0; newTimes = find(binDum ==1 ); cnt=cnt+1 ;
%            figure;hold on; plot(binary,'k+'); plot(binDum,'r+');plot(hit,1.1,'b*');
%            xlim([rebound(1)-15 rebound(end)+5]);ylim([0.96 1.12]);pause; clf
        end
% % Bin "donwstream neuron" spiketrain at larger binsize (to match regular convergent model binning)
    spikeTrain = find( binDum ==1) .*binSize_sm ; %!!! modeled spikeTrain (converted to spike times in sec)
    dum = histc(spikeTrain,rasterTime); boop(tr,:) = dum(1:end-1) ; % BINNED RASTER
    
    spikeTrains{s}{f}{tr} = spikeTrain ; 
%figure;hold on ; plot(merged(tr,:),'k');plot(epspDum,'r','linewidth',2);plot(find(binary==1),3,'b.','markersize',8)
    end % END (tr) per trial
    spykes{s}{f} = boop; % EPSP+threshold model, binned spikerate {stimulus} (trial x time) 
    templates{s}(f,:) = mean(boop); % means for each fMod
    
    end % END (s) per shuffle
end % END (f) per modulation frequency [[raster construction]]
    
% % % Decode stimulus per shuffle % % %     
for ss = 1:nShuffles
conf = zeros(nFMs,nFMs); % initialize confusion matrix
for fM = 1:nFMs % for each mod freq
    hits = zeros(1,nFMs); % intialize single row of confusion matrix
for i = 1:size(spykes{ss}{fM},1) % for each trial
    dum = spykes{ss}{fM};  % all trials for (fM)
    dum(i,:)=[]; % remove (i) trial for Xvalidation
    temp_xV = templates{ss}; temp_xV(fM,:) = mean(dum) ; % update templates with (i) trial exclusion
    eucDists = dist(spykes{ss}{fM}(i,:),temp_xV') ; % dist btwn "tr" trial and all others
    decoded = find(eucDists == min(eucDists)) ; 
% % % % WHAT TO DO WHEN MORE THAN ONE STIM IS AT MIN??? (i.e. equal prob.) % % % %
    hits(decoded) = hits(decoded) + 1/length(decoded); % IFF more than one, split the difference
end % end (i) per trial confustion matrixs
conf(fM,:) = hits ./ i ; % normalize to yield probability

end % % end (fM)
conf_all{ss} = conf ;
performance(ss) = sum(diag(conf)) / length(diag(conf)) ;
end % END (ss) per shuffle

end
