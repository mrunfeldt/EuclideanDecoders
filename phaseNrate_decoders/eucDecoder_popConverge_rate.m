% % Rate Only, Population Convergent Model, Euclidean Distance decoder % %

% % % Rate Only, Convergent Model: For each trial, sums neuronal activity of
% all neurons over entire duration of stimulus % % %

% % For rate classifier, Convergent Model is the same as the Labled Line
% Model % % :-O

% INPUT (1) : SPIKETIMES (in sec):  data{neuron}{stimulus}{trial}

% % MJRunfeldt 2015_05_28

function[performance,conf,rasters] = eucDecoder_popConverge_rate(input)

%input = spykes(us); input=data;% {ch}{fM}{tr}


nFMs = mode(cellfun(@length,input)) ; % number of modulation frequencies
% % % Bin spike times and get all-trial template means % % 
rasters = cell(1,nFMs) ;templates = [] ;
for f = 1:nFMs % for each mod freq
% % % Sum NEURONS IN TIME, "Convergence" % % %    
    % requires channels to have same duration and # of trials%
    nTrials = length(input{1}{f}); 
        merged = zeros(1,nTrials);
    for nn = 1:length(input) % for each neuron
    	raz = cellfun(@length,input{nn}{f}) ; % # of spikes per trial
        merged = merged+raz; % nTrials x Time
    end
    rasters{f} = merged ; 
    templates(f,:) = mean(merged); % means for each ModF
end


conf = zeros(nFMs,nFMs); % initialize confusion matrix
for fM = 1:nFMs % for each mod freq
    hits = zeros(1,nFMs); % intialize single row of confusion matrix
for i = 1:size(rasters{fM},1) % for each trial
    dum = rasters{fM};  % all trials for (fM)
    dum(i)=[]; % remove (i) trial for Xvalidation
    temp_xV = templates; temp_xV(fM,:) = mean(dum) ; % update templates with (i) trial exclusion
    eucDists = dist(rasters{fM}(i),temp_xV') ; % dist btwn "tr" trial and all others
    decoded = find(eucDists == min(eucDists)) ; 
% % % % WHAT TO DO WHEN MORE THAN ONE STIM IS AT MIN??? (i.e. equal prob.) % % % %
    hits(decoded) = hits(decoded) + 1/length(decoded); % IFF more than one, split the difference
end % end (i) per trial confustion matrix
conf(fM,:) = hits ./ i ; % normalize to yield probability
%figure;imagesc(conf)
end % % end (fM)
performance = sum(diag(conf)) / length(diag(conf)) ;

end