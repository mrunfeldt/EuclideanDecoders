% % % FULL SPIKE TRAIN CLASSIFER: FOR NEURONAL PAIRS % % %
% Euclidean distance decoder of modulation frequency: Convergent Model population decoder % % %
% % % Convergent Model == sums activity of neurons over time % % %
% % % % Removes "noise" trial-trial correlations by reshuffling neurons
% across trials while maintaining perstimulus timing % % % %
% % % % Also records pairwise correlations in binned spike count for total and
% with trial-trial corr's removed: DESIGNED EXPLICITLY FOR TWO NEURONS
% % % % Pairwise Correlations obtained per FM, using binned counts
% vectorized over trials % % % % % % % % % % % % % % % % % % % % % % % % %

% INPUT (1) : SPIKETIMES (in sec):  data{nrn}{fm}{tr} 
% (2): binSize for rasters (same sampling rate as spiketimes)
% (3): number of shuffles for decorrelation
% (4): duration of stimulus (in sec)
% % MJRunfeldt 05_11_2015

function[performance,conf_all,pairCorr,rMTF,tMTF] = eucDecoder_pairsConverge_dcShuffle_wCorrs_wMTF(input,binSize,nShuffles,duration)

%input=data;% nShuffles = 5; input=data;input = spykes(us); % {ch}{fM}{tr} % for troubleshooting

time = 0:binSize:duration;
nFMs = mode(cellfun(@length,input)) ; % number of modulation frequencies
nChan = length(input) ; % # of channels
nTrials = length(input{1}{1}) ; %  # of trials (assumes nTrials is equal for all neurons)  

% % % % % % % % % FOR each shuffled population, calculate MTF and generate rasters
% for decoding % % % % % % % % %

pairCorr = zeros(nFMs,nShuffles) ; performance = zeros(1,nShuffles); conf_all=cell(1,nShuffles);
for s = 1:nShuffles % for each shuffling

    convRaster = cell(1,nFMs) ; % convergent rasters {fm}(tr x time)
    convSpikes = cell(1,nFMs) ; % convergent spike trains {fm}{tr}
    templates = zeros(nFMs,length(time)-1) ; % convergent raster trial-mean (fm,time)

% % % Randomly permute trials across modulation frequences and neurons (independetly)        
newNrn = cell(1,nChan) ; % {nn}{fm}{tr}
    for nn = 1:nChan % for each channel/neuron,randomly permute trials
        newNrn{nn} = cell(1,nFMs) ; % {chan}{fm}{tr}
        for fm = 1:nFMs % for each modulation freq
            newNrn{nn}{fm} = input{nn}{fm}(randperm(nTrials)) ; % {nn}{fm}{tr}
        end % END (fm) mod freq
    end % END (nn) neuron
        
% Generate "convSpikes" spiketimes for MTF and "convRasters" rasters for decoding %   
% % Also obtain pairwise Corr's for each FM (vectorize binned counts over trials)
    for fm = 1:nFMs % for each modulation freq
        convSpikes{fm} = cell(1,nTrials) ;  
% % % for each fm/trial, generate convergent spiketimes % % %
        for tr = 1:nTrials
            for nn = 1:nChan
                convSpikes{fm}{tr} = [convSpikes{fm}{tr} newNrn{nn}{fm}{tr}];
            end
        convSpikes{fm}{tr} = sort(convSpikes{fm}{tr},'ascend'); % resort (even though it doesn't matter if you are binning)
            dum = histc(convSpikes{fm}{tr},time) ; % bin spike counts in time 
            if isempty(dum); dum = zeros(1,length(time));end % Add zeros if no spikes
        convRasters{fm}(tr,:) = dum(1:length(time)-1) ; % remove that pesky zero bin
        end % END (tr) 
        templates(fm,:) = mean(convRasters{fm}) ; % average over trials (length time)
            
% % % CALCULATE PAIRWISE CORRELATIONS - 1 per FM - vectorize over trials %
    nrnBin = cell(1,nChan) ; 
        for nn = 1:nChan
            binnedNrn = [];
            for tr = 1:nTrials
                mySpikes = histc(newNrn{nn}{fm}{tr},time); % bin spikes
                if isempty(mySpikes);mySpikes=zeros(1,length(time));end % add zeros if no spikes
                binnedNrn = [binnedNrn mySpikes(1:end-1)]; % concat trials
            end
            nrnBin{nn} = binnedNrn ; 
        end % END (nn) 
        [rDum]= corrcoef(nrnBin{1},nrnBin{2}) ; 
        pairCorr(fm,s) = rDum(2) ; 
    end % END (fm) per FM (within obtain rasters and corrs loop)
    

% % % % % DECODE MODULATION FREQUENCY % % % % % % %
conf = zeros(nFMs,nFMs); % initialize confusion matrix
for fM = 1:nFMs % for each mod freq
    hits = zeros(1,nFMs); % intialize single row of confusion matrix
for i = 1:nTrials % for each trial
    dum = convRasters{fM};  % all trials for (fM)
    dum(i,:)=[]; % remove (i) trial for Xvalidation
    temp_xV = templates; temp_xV(fM,:) = mean(dum) ; % update templates with (i) trial exclusion
    eucDists = dist(convRasters{fM}(i,:),temp_xV') ; % dist btwn "tr" trial and all others
    decoded = find(eucDists == min(eucDists)) ; 
% % % % WHAT TO DO WHEN MORE THAN ONE STIM IS AT MIN??? (i.e. equal prob.) % % % %
    hits(decoded) = hits(decoded) + 1/length(decoded); % IFF more than one, split the difference
end % end (i) per trial confustion matrix
conf(fM,:) = hits ./ i ; % normalize to yield probability
end % % end (fM)
% figure;imagesc(conf);colorbar;xlabel('Decoded');ylabel('Actual')

performance(s) = sum(diag(conf)) / length(diag(conf)) ;
conf_all{s} = conf; 

end % END (s)

% % % % % % MTFs are the SAME no matter how you shuffle (as long as the
% cumulative trials are the same % % %  % % % % % % % % % % % % % % % % % 
% % % Calculate MTF % % % % input for MTF is: {fM}{trial}    
[rMTF,tMTF] = create_MTF(convSpikes,[1:nFMs],duration) ; % MTF for each shuffled pop    

end % END main function 


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % MTF subroutine function % % % 
function [rMTF,tMTF] = create_MTF(raw,fMs,duration,varargin)

% raw = convSpikes ; fMs = sorted.fMs ; raw = singleChan; % % for troublshoot
%disp('yo - lin 110')
for f = 1:length(fMs)
    modfreq = fMs(f) ; % in secnds
    allSpikes = cell2mat(raw{f}) ; % all spike times all trials
    tMTF(f) = sqrt(sum(sin(2*pi*modfreq*allSpikes))^2 + ...
            sum(cos(2*pi*modfreq*allSpikes))^2)/length(allSpikes);
    
    nTrials = length(raw{f});
    rMTF(f) = length(allSpikes) / (nTrials*duration) ;

end

if nargin > 3 & varargin{1} == 1 % plot only tMTF
f1=figure;plot(log10(fMs),tMTF,'.-k','linewidth',3)
set(gca,'xtick',log10(fMs),'xticklabel',num2cell(fMs))
xlabel('Modulation Frequency (Hz)');ylabel('Vector Strength')
title('Temporal Modulation Transfer Function');shg

elseif nargin > 3 & varargin{1} == 2 % plot tMTF and rMTF
f1=figure;plot(log10(fMs),tMTF,'.-k','linewidth',3)
set(gca,'xtick',log10(fMs),'xticklabel',num2cell(fMs))
xlabel('Modulation Frequency (Hz)');ylabel('Vector Strength')
title('Temporal Modulation Transfer Function');shg    
    
f2=figure;plot(log10(fMs),rMTF,'.-k','linewidth',3)
set(gca,'xtick',log10(fMs),'xticklabel',num2cell(fMs))
xlabel('Modulation Frequency (Hz)');ylabel('Firing Rate (Hz)')
title('Rate Modulation Transfer Function');shg
end
end