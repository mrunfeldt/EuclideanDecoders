
% INUT: (1) Raw Spike Times: either cell for each trial or single trial
% (2) binSize in whatever sampling input is in (e.g. spiketimes in ms, binsize in ms)
% (3) Optional, duration of trial

% % % bins rasters into vector (single trial) or array (multiple trials)
% MJRunfeldt 2015_02_09

% stIn = input{1}; 
% binSize = 5e-3 ;% binSize in whatever sampling input is in (e.g. spiketimes in ms, binsize in ms)
% %trialDur = 1 ; 

function [rasters,time] = spikeTimes_toRasters(stIn,binSize,varargin)
% stIn = input{f}; trialDur=duration;

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

end