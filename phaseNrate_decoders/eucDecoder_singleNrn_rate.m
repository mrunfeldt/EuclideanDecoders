% % % "RATE-ONLY" CLASSIFYER: response is the sum of activity across entire
% duration of stimulus % % %

% % Euclidean distance decoder for single neuron or channel: decodes each
% trial based on vector of firing rate in binned time. Removes decoded
% trial from stimulus mean for cross-validation %

% INPUT (1) : SPIKETIMES (in sec):  data{channel}{stimulus}{trial} 
% (2): binSize for rasters (same sampling rate as spiketimes)
% (3) duration of stimulus
% % MJRunfeldt 2015_05_27

function[performance,conf,FR] = eucDecoder_singleNrn_rate(input)

% input = singleChan; input = SU;
nFMs = length(input) ; % number of modulation frequencies (stimuli)
nTRs = length(input{1}) ; % # number of trials

% % % Bin spike times and get all-trial template means % % 
FR = zeros(nFMs,nTRs) ;templates = [] ;
for f = 1:nFMs % for each mod freq
    FR(f,:) = cellfun(@length,input{f}) ; % N = # of trials
    templates(f) = mean(FR(f,:)) ; % mean rate for each ModF
end

conf = zeros(nFMs,nFMs); % initialize confusion matrix
for fM = 1:nFMs % for each mod freq
    hits = zeros(nFMs,1); % intialize single row of confusion matrix
for i = 1:nTRs % for each trial
    dum = FR(fM,:);  % all trials for (fM)
    dum(i)=[]; % remove (i) trial for Xvalidation
    temp_xV = templates;  % update templates with (i) trial exclusion
    temp_xV(fM) = mean(dum)  ; % 
    eucDists = dist(FR(fM,i),temp_xV) ; % dist btwn "tr" trial and each modulation Freq
    decoded = find(eucDists == min(eucDists)) ; 
% % % % WHAT TO DO WHEN MORE THAN ONE STIM IS AT MIN??? (i.e. equal prob.) % % % %
    hits(decoded) = hits(decoded) + 1/length(decoded); % IFF more than one, split the difference
end % end (i) per trial confustion matrix
conf(fM,:) = hits ./ i ; % normalize to yield probability
end % % end (fM)
%figure;imagesc(conf);xlabel('decoded');ylabel('actual')

performance = sum(diag(conf)) / length(diag(conf)) ;

end