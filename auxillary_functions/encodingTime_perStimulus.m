
% % Calculates Encoding Time For each stimulus by (1) taking the autocorrelation of a
% neuron's psth (2) fitting a gaussian function via nonlinear least squares fitting
% (3) EncodingTime = 2 * std(gausFit)

% INPUT FORMAT: (1) "input": spikeTimes in seconds {stimulus}{trial}
% (2) duration of stimulus (in seconds)
% OUTPUT: Encoding time in ms (1 per stimulus)

% MJRunfeldt 2015_07_06


function [ET]  = encodingTime_perStimulus(input,duration)
% input = spikeTimes ; % For Devel

% % HARDCODED parameters % % 
binSize = 1e-3 ; %  1 ms
window = 30e-3; %  30 ms : max lag of autocorr
% Gaussian function: parameters - (1) height (2) mean (3) std
gausFun = @(hms,x) hms(1) .* exp (-(x-hms(2)).^2 ./ (2*hms(3)^2)) ; 
% % % % % % % % % % % % % % %


% % Obtain psth and EncodingTime per stimulus
psth = []; 
for nFm = 1:length(input)  % rasters for population coupling
    
    if length(input{nFm}) < 5 ; ET(nFm) = 0 ; % set to zero if not at least X spikes
    else
        dum = spikeTimes_toRasters(input{nFm},binSize,duration) ;% trial x Time
        psth = sum(dum,1) ; 

        [cross,lagz]=xcorr(psth,psth,round(window/binSize));
        zToss = find(lagz==0); cross(zToss)=[];lagz(zToss)=[]; % REMOVE ZERO

        cross = mat2gray(cross) ; % normalize auto-correlation
        
x0 = [1; 0; round(window/binSize)/2 ]; % initial estimates of paramters

options=optimset('Display','off'); % turn off lsq messages
[params]=lsqcurvefit(gausFun,x0,[1:length(lagz)],cross,[],[],options); % meat and potatoes
sigma = (params(3)/binSize)*1e-3 ; % in ms

        ET(nFm) =  2*sigma   ; % encoding time in ms
    end
     % figure;plot(lagz,cross);title(num2str(ET(nFm))); pause;close
end
    

end % END main function



% % % Subroutine % % %
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
time(end) = []; 

end
