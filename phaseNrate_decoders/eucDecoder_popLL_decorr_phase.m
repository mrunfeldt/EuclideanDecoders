% % % PHASE - ONLY CLASSIFYIER: Labeled Line Model, decorrelated, Euclidean distance decoder

% % % Labeled Line,  Phase Only : Concats neurons in time to create a single vector
% length [nTimeBins * nNeurons], then normalizes by mean spike rate over
% entire vector to remove rate information. Does NOT normalize per neuron.

% % DECORRELATED: Reshuffles across trials to remove "noise" correlations % % 
% % % % OUTPUT includes decoding performance, confusion matricies, and
% rasters for each trial shuffling % % % 

% INPUT (1) : SPIKETIMES (in sec):  data{neuron}{stimulus}{trial} 
% (2): binSize for rasters (same sampling rate as spiketimes)
% % MJRunfeldt 2015_05_28

function[performance,conf_all,rasters] = eucDecoder_popLL_decorr_phase(input,binSize,nShuffles,duration)

%nShuffles = 5; input = data ;
%input = dataIn(nrns); % {ch}{fM}{tr}
time = 0:binSize:duration;

nFMs = mode(cellfun(@length,input)) ; % number of modulation frequencies
nChan = length(input) ; % # of channels
% % % Bin spike times and get all-trial template means % % 
rasters = cell(1,nFMs) ;templates = [] ;
for f = 1:nFMs % for each mod freq
% % % CONCAT NEURONS IN TIME "labeled line" % % %    
    % requires channels to have same duration and # of trials%
     merged = [];
    for nn = 1:nChan % for each channel/neuron
        [raz] = spikeTimes_toRasters(input{nn}{f},binSize,duration) ; % Bin Spike Times
        merged = horzcat(merged,raz); % nTrials x Time
    end

% % % Generate multiple shuffles of channel-trial simultaneity % % %    
% % Phase - only : normalize after shuffling
ori = merged ; nTrials = size(ori,1) ; % original: [trials X catTime]
oneTrial = size(raz,2) ; % time in single trial per channel
for s = 1:nShuffles    
     timeTrack = 1; new =[];
    for ch = 1:nChan
        singleC = ori(:,timeTrack:timeTrack + oneTrial - 1) ;
        mix = singleC(randperm(nTrials),:);  
        %figure;imagesc( vertcat(mean(singleC),mean(mix)));colorbar
        new = horzcat(new,mix) ;
        timeTrack = timeTrack + oneTrial ; 
    end ;% figure;imagesc( vertcat(mean(new), mean(merged))); title('Trial means are same')
    dum=bsxfun(@rdivide,new',mean(new'))'; dum(isnan(dum)) = 0; % NORMALIZE by mean, remove NaNs
    rasters{s}{f} = dum;
    templates{s}(f,:) = mean(dum); % means for each fM
end % end (s) shuffles
end % end (f) 

%figure;imagesc(cell2mat(cellfun(@mean,templates,'un',0)'));title('Collapsed templates same for all shuffles')

% % % Generate CONFUSION MATRIX per trial- shuffling % % % 
%conf_mean=zeros(nFMs,nFMs); % initialize confusion matrix
conf_all = cell(1,nShuffles) ;
for s = 1:nShuffles

conf = zeros(nFMs,nFMs); % initialize confusion matrix
for fM = 1:nFMs % for each mod freq
    hits = zeros(1,nFMs); % intialize single row of confusion matrix
for i = 1:size(rasters{s}{fM},1) % for each trial
    dum = rasters{s}{fM};  % all trials for (fM)
    dum(i,:)=[]; % remove (i) trial for Xvalidation
    temp_xV = templates{s}; temp_xV(fM,:) = mean(dum) ; % update templates with (i) trial exclusion
    eucDists = dist(rasters{s}{fM}(i,:),temp_xV') ; % dist btwn "tr" trial and all others
    decoded = find(eucDists == min(eucDists)) ; 
% % % % WHAT TO DO WHEN MORE THAN ONE STIM IS AT MIN??? (i.e. equal prob.) % % % %
    hits(decoded) = hits(decoded) + 1/length(decoded); % IFF more than one, split the difference
end % end (i) per trial confustion matrix
conf(fM,:) = hits ./ i ; % normalize to yield probability
%figure;imagesc(conf);xlabel('Decoded');ylabel('Actual')
end % % end (fM)

conf_all{s} = conf; 
%conf_mean= conf_mean+conf;
performance(s) = sum(diag(conf)) / length(diag(conf)) ;
end % end (s) shuffle
%conf_mean= conf_mean./s;
end