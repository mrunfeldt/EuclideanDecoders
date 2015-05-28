% % Full SpikeTrain Classifier for single neuron % % 

% % Euclidean distance decoder for single neuron or channel: decodes each
% trial based on vector of firing rate in binned time. Removes decoded
% trial from stimulus mean for cross-validation %

% INPUT (1) : SPIKETIMES {neuron}{stimulus}{trial} 
% (2): binSize for rasters (same sampling rate as spiketimes)
% (3) duration of stimulus (same sampling rate as spiketimes)
% % MJRunfeldt 02_09_2015

function[performance,conf,rasters] = eucDecoder_singleNrn_full(input,binSize,duration)

% input = singleChan; input = MU;
nFMs = length(input) ; % number of modulation frequencies
time = 0:binSize:duration;

% % % Bin spike times and get all-trial template means % % 
rasters = cell(1,nFMs) ;templates = [] ;
for f = 1:nFMs % for each mod freq
    [rasters{f}] = spikeTimes_toRasters(input{f},binSize,duration) ; % Bin Spike Times
    templates(f,:) = mean(rasters{f}); % means for each ModF
end

conf = zeros(nFMs,nFMs); % initialize confusion matrix
for fM = 1:nFMs % for each mod freq
    hits = zeros(1,nFMs); % intialize single row of confusion matrix
for i = 1:size(rasters{fM},1) % for each trial
    dum = rasters{fM};  % all trials for (fM)
    dum(i,:)=[]; % remove (i) trial for Xvalidation
    temp_xV = templates; temp_xV(fM,:) = mean(dum) ; % update templates with (i) trial exclusion
    eucDists = dist(rasters{fM}(i,:),temp_xV') ; % dist btwn "tr" trial and all others
    decoded = find(eucDists == min(eucDists)) ; 
% % % % WHAT TO DO WHEN MORE THAN ONE STIM IS AT MIN??? (i.e. equal prob.) % % % %
    hits(decoded) = hits(decoded) + 1/length(decoded); % IFF more than one, split the difference
end % end (i) per trial confustion matrix
conf(fM,:) = hits ./ i ; % normalize to yield probability
end % % end (fM)
%figure;imagesc(conf)

performance = sum(diag(conf)) / length(diag(conf)) ;

end