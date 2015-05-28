% % % PHASE-ONLY single neuron classifier: Normalizes each trial and template means so
% that only temporal information is used to decode stimulus % % %

% % Euclidean distance decoder for single neuron or channel. Removes decoded
% trial from stimulus mean for cross-validation %

% INPUT (1) : SPIKETIMES (in sec):  data{channel}{stimulus}{trial} 
% (2): binSize for rasters (same sampling rate as spiketimes)
% (3) duration of stimulus
% % MJRunfeldt 2015_05_27

function[performance,conf,rasters] = eucDecoder_singleNrn_phase(input,binSize,duration)

% input = singleChan; input = SU;
nFMs = length(input) ; % number of modulation frequencies (stimuli)
time = 0:binSize:duration;

% % % Bin spike times and get all-trial template means % % 
rasters = cell(1,nFMs) ;templates = [] ;
for f = 1:nFMs % for each mod freq
    [raz] = spikeTimes_toRasters(input{f},binSize,duration) ; % Bin Spike Times
    dum=bsxfun(@rdivide,raz',mean(raz')); dum(isnan(dum)) = 0; % Normalize by mean, remove NaNs
    rasters{f} = dum' ; % trial x time
    templates(f,:) = mean(rasters{f}) ./ mean(mean(rasters{f})); % means for each ModF
end

conf = zeros(nFMs,nFMs); % initialize confusion matrix
for fM = 1:nFMs % for each mod freq
    hits = zeros(1,nFMs); % intialize single row of confusion matrix
for i = 1:size(rasters{fM},1) % for each trial
    dum = rasters{fM};  % all trials for (fM)
    dum(i,:)=[]; % remove (i) trial for Xvalidation
    temp_xV = templates;  % update templates with (i) trial exclusion
    temp_xV(fM,:) = mean(dum) ./ mean(mean(dum)) ; % don't forget to NORMALIZE (remove rate info)
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