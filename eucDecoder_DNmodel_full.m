% % % Decoding Performance of "Downstream Neuron" model % % % 
% Bin spike times at high temporal resolution and then convolves binary signal
% with decaying exponential to simulate neural EPSPs.
% Designed for multineuronal Convergent model, or "Pooled" Model.
% % Applies threshold to EPSPs to generate downstream neuron model; includes (hardcoded) refractory period.
% % Decoding is based on binned raster generated from EPSP model

% INPUT: (1) spikeTimes {neuron}{stimulus}{trial}
% (2) tau (~full width of EPSP @ half height)
% (3) ISI rebound (4) EPSP threshold (will be converted to fnc of EPSP amp)
% (4) small binSize for EPSP
% (5) binSize for output raster (6) duration of stimulus 
% % Time unit for all of input must be same (just use seconds)

% OUTPUT: (1) decoding performance (2) confusion matrix 
% (3) binned raster of DN model (firing rates) - {stimulus} (trial x time) 
% (4) spike times of DN model (time of each spike) :  {stimulus} {trial} (spikeTime) 

% MJRunfeldt 2015_08_14: Updated EPSP model to use single eponential,
% Threshold defined relative to maximum amplitude of EPSP, which is:
% Tau / 6.6589

function [performance, conf, spykes, spikeTrains] = ...
    eucDecoder_DNmodel_full(data,tau,rebound,threshold,binSize_sm,binSize_ras,duration)

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

% % % Bin spike times and get all-trial template means % % 
nFMs = mode(cellfun(@length,data)) ; % number of modulation frequencies
rasters = cell(1,nFMs) ;spykes = rasters; templates = [] ; spikeTrains=[];
%EPSPs = rasters;
for f = 1:nFMs % for each mod freq
% % % Sum NEURONS IN TIME, "Convergence" % % %    
    % requires channels to have same duration and # of trials%
    nTrials = length(data{1}{f}); 
        merged = zeros(nTrials,length(time)-1);
    for nn = 1:length(data) % for each neuron
    	raz = spikeTimes_toRasters(data{nn}{f},binSize_sm,duration) ; % Bin Spike Times SECONDS
        if isempty(raz); raz = zeros(1,length(time)-1);end
        merged = merged+raz; % nTrials x Time
        
    end % END (nn)
    rasters{f} = merged; % binned rasters {fm} (trial x time) 
    
    boop=[]; % epsp vector (i.e. binary signal convolved with decaying exponential)
    for tr=1:nTrials       
        epspDum = conv(merged(tr,:),myExp); epspDum = epspDum(1:size(merged,2));
        %EPSPs{f}{tr} = epspDum; 
        binary = double(epspDum>thresh) ;% binary spiketrain
        spikeTimes = find(binary==1); binDum = binary ; newTimes = spikeTimes; cnt=1;
% % march forward through binary and set "isiLim" to 0 after each spike
        while cnt < length(newTimes) 
           hit = newTimes(cnt) ; rebound = hit+1 : hit+1+isiLim ; % time frame of rebound
           binDum(rebound)=0; newTimes = find(binDum ==1 ); cnt=cnt+1 ;
%            figure;hold on; plot(binary,'k+'); plot(binDum,'r+');plot(hit,1.1,'b*');
%            xlim([rebound(1)-15 rebound(end)+5]);ylim([0.96 1.12]);pause; clf
        end
    spikeTrain = find( binDum ==1) .*binSize_sm ; %!!! modeled spikeTrain (converted to spike times in sec)
    dum = histc(spikeTrain,rasterTime); boop(tr,:) = dum(1:end-1) ; % BINNED RASTER
    
    spikeTrains{f}{tr} = spikeTrain;
%figure;hold on ; plot(merged(tr,:),'k');plot(epspDum,'r','linewidth',2);
%plot(find(binary==1),3,'b.','markersize',6); plot(find(binDum==1),3,'g.','markersize',13);
    end % END (tr) per trial
    spykes{f} = boop; % binned spiketrain from EPSP+thresh model {stimulus} (trial x time) 
    templates(f,:) = mean(boop); % means for each fMod
end
    
    
conf = zeros(nFMs,nFMs); % initialize confusion matrix
for fM = 1:nFMs % for each mod freq
    hits = zeros(1,nFMs); % intialize single row of confusion matrix
for i = 1:size(spykes{fM},1) % for each trial
    dum = spykes{fM};  % all trials for (fM)
    dum(i,:)=[]; % remove (i) trial for Xvalidation
    temp_xV = templates; temp_xV(fM,:) = mean(dum) ; % update templates with (i) trial exclusion
    eucDists = dist(spykes{fM}(i,:),temp_xV') ; % dist btwn "tr" trial and all others
    decoded = find(eucDists == min(eucDists)) ; 
% % % % WHAT TO DO WHEN MORE THAN ONE STIM IS AT MIN??? (i.e. equal prob.) % % % %
    hits(decoded) = hits(decoded) + 1/length(decoded); % IFF more than one, split the difference
end % end (i) per trial confustion matrix
conf(fM,:) = hits ./ i ; % normalize to yield probability
%figure;imagesc(conf)
end % % end (fM)
performance = (sum(diag(conf)) / length(diag(conf))) *100 ; % percent correct

end % END main function



% % % % Subroutine % % %
function [rasters,time] = spikeTimes_toRasters(stIn,binSize,varargin)
% stIn = data; trialDur=duration;

% % % Assign trial duration if variable not defined: dur = ceil(last spike time)
if nargin == 2 & iscell(stIn)% if no trial duration stated, extract
    trialDur = ceil( max( cell2mat( cellfun(@(x) x(end),stIn,'un',0) ) ) ) ; % length of each
elseif nargin ==2 
    trialDur = ceil( stIn(end) ) ;
elseif nargin == 3
    trialDur = varargin{1} ;
end

time = 0: binSize : trialDur ; % % BIN TIME

if iscell(stIn) % if input is cell, assumes each cell is a trial
    rasters = zeros(length(stIn),length(time)-1) ; % trials x time
for tr = 1:length(stIn) % for each trial
    if ~isempty(stIn{tr})
    dum = histc(stIn{tr},time);dum(end)=[]; % removes last bin because never filled
    rasters(tr,:) = dum ;
    else rasters(tr,:) = zeros(1,length(time)-1) ;
    end
end
else % otherwise input is single trial array
    dum = histc(stIn,time);dum(end)=[]; % removes last bin because never filled
    rasters = dum ;
end
time(end) = []; 

end
