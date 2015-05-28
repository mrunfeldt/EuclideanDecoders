% % Full Spike Train, Population Convergent Model euclidean decoder % % 

% % % Convergent Model: Sums neuronal activity of population over time. 

% INPUT (1) SPIKETIMES (in sec):  data{neuron}{stimulus}{trial}
% (2): binSize for rasters (same sampling rate as spiketimes)
% (3) Duration of stimulus 
% % MJRunfeldt 02_09_2015

function[performance,conf,rasters] = eucDecoder_popConverge_full(input,binSize,duration)

%input = spykes(us); % {ch}{fM}{tr} % for Devel

time = 0:binSize:duration;

nFMs = mode(cellfun(@length,input)) ; % number of modulation frequencies
% % % Bin spike times and get all-trial template means % % 
rasters = cell(1,nFMs) ;templates = [] ;
for f = 1:nFMs % for each mod freq
% % % Sum NEURONS IN TIME, "Convergence" % % %    
    % requires channels to have same duration and # of trials%
    nTrials = length(input{1}{f}); 
        merged = zeros(nTrials,length(time)-1);
    for nn = 1:length(input) % for each neuron
    	raz = spikeTimes_toRasters(input{nn}{f},binSize,duration) ; % Bin Spike Times
        merged = merged+raz; % nTrials x Time
    end
    rasters{f} = merged; 
    templates(f,:) = mean(merged); % means for each ModF
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
%figure;imagesc(conf)
end % % end (fM)
performance = sum(diag(conf)) / length(diag(conf)) ;

end