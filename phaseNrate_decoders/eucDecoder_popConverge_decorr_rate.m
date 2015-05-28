% % % RATE - ONLY CLASSIFYIER: Convergent Model, decorrelated, euclidean
% distance decoder % % % 

% % % Rate Only Convergent Model: For each trial, sums neuronal activity of
% all neurons over entire duration of stimulus % % %

% % For rate classifier, Convergent Model is the same as the Labled Line
% Model % % :-O

% % DECORRELATED: Reshuffles across trials to remove "noise" correlations % % 
% % % % OUTPUT includes decoding performance, confusion matricies, and
% rasters for each trial shuffling % % % 

% INPUT (1) : SPIKETIMES (in sec):  data{neuron}{stimulus}{trial}
% (2) : number of shuffles for decorrelation
% % MJRunfeldt 2015_05_27

function[performance,conf_all,rasters] = eucDecoder_popConverge_decorr_rate(input,nShuffles)

% nShuffles = 5; % input = spykes(us); % {ch}{fM}{tr} % For devel

nFMs = mode(cellfun(@length,input)) ; % number of modulation frequencies
nChan = length(input) ; % # of channels
% % % Bin spike times and get all-trial template means % % 
rasters = cell(1,nChan) ;templates = [] ;
for f = 1:nFMs % for each mod freq
    
% % % Get Rasters in binned Time % % %    
    for nn = 1:nChan % for each channel/neuron
        rasters{nn} = cellfun(@length,input{nn}{f}); % # of spikes. Vector nTrials long
    end

% % % Generate CONVERGENT model from multiple shuffles of channel-trial simultaneity % % %    
nTrials = length(input{1}{f}) ; 
for s = 1:nShuffles    
     merged = zeros(size(rasters{1})) ;
     for ch = 1:nChan
         singleC = rasters{ch} ; % original raster
         merged = merged + singleC(randperm(nTrials));  % mix trials
     end % figure;imagesc( vertcat(mean(new), mean(merged))); title('Trial means are same')
    shufflRasters{s}{f} = merged; % one per trial
    templates{s}(f,:) = mean(merged); % means for each ModF
end % end (s) shuffles
end % end (f) 

%figure;imagesc(cell2mat(cellfun(@mean,templates,'un',0)'));title('Collapsed templates same for all shuffles')

% % % Generate CONFUSION MATRIX per trial- shuffling % % % 
conf_all = cell(1,nShuffles) ; 
%conf_mean=zeros(nFMs,nFMs); % initialize confusion matrix
for s = 1:nShuffles

conf = zeros(nFMs,nFMs); % initialize confusion matrix
for fM = 1:nFMs % for each mod freq
    hits = zeros(1,nFMs); % intialize single row of confusion matrix
for i = 1:size(shufflRasters{s}{fM},1) % for each trial
    dum = shufflRasters{s}{fM};  % all trials for (fM)
    dum(i)=[]; % remove (i) trial for Xvalidation
    temp_xV = templates{s}; temp_xV(fM,:) = mean(dum) ; % update templates with (i) trial exclusion
    eucDists = dist(shufflRasters{s}{fM}(i),temp_xV') ; % dist btwn "tr" trial and all others
    decoded = find(eucDists == min(eucDists)) ; 
% % % % WHAT TO DO WHEN MORE THAN ONE STIM IS AT MIN??? (i.e. equal prob.) % % % %
    hits(decoded) = hits(decoded) + 1/length(decoded); % IFF more than one, split the difference
end % end (i) per trial confustion matrix
conf(fM,:) = hits ./ i ; % normalize to yield probability
end % % end (fM)
%figure;imagesc(conf);colorbar;xlabel('Decoded');ylabel('Actual')

performance(s) = sum(diag(conf)) / length(diag(conf)) ;
conf_all{s} = conf; 
%conf_mean = conf_mean +conf;
end % end (s) shuffle
%conf_mean = conf_mean ./ s ; % Mean of all confusion matricies

end