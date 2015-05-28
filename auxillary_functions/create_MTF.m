
% % % Creates temporal and rate Modulation Transfer Function for a single neuron/channel % % % 
% tMTF is vector strength vs modulation frequency. 
% rMTF is average number of spikes (per trial) vs modulation frequency. 
% INPUTS:(1) "raw" : spikeTimes in seconds {fM}{trial}
% (2) "fMs" : vector of modulation frequencies (in Hz)
% (3) "duration" : duration of stimulus in seconds
% (4) varargin: set to 1 if you want to plot tMTF, set to 2 to plot both
% tMTF and rMTF. Else leave empty
% INPUT IS ASSUMED IN SECONDS

% Created MJRunfeldt April 2015

function [rMTF,tMTF] = create_MTF(raw,fMs,duration,varargin)

% raw = newSpikes ; fMs = sorted.fMs ; raw = singleChan; % % for troublshoot


for f = 1:length(fMs)
    modfreq = fMs(f) ; % in secnds
    
    try
        allSpikes = cell2mat(raw{f}) ; % all spike times all trials
    catch
        allSpikes = cell2mat(raw{f}') ; % cell2mat only likes [1 x N] not [N x 1]
    end
        
    tMTF(f) = sqrt(sum(sin(2*pi*modfreq*allSpikes))^2 + ...
            sum(cos(2*pi*modfreq*allSpikes))^2)/length(allSpikes);
    
    nTrials = length(raw{f});
    rMTF(f) = length(allSpikes) / (nTrials*duration) ;

end

if nargin > 3 & varargin{1} == 1 % plot only tMTF
f1=figure;plot(log10(fMs),tMTF,'.-k','linewidth',3)
set(gca,'xtick',log10(fMs),'xticklabel',num2cell(fMs))
xlabel('Modulation Frequency (Hz)');ylabel('Vector Strength')
title('Temporal Modulation Transfer Function');shg

elseif nargin > 3 & varargin{1} == 2 % plot tMTF and rMTF
f1=figure;plot(log10(fMs),tMTF,'.-k','linewidth',3)
set(gca,'xtick',log10(fMs),'xticklabel',num2cell(fMs))
xlabel('Modulation Frequency (Hz)');ylabel('Vector Strength')
title('Temporal Modulation Transfer Function');shg    
    
f2=figure;plot(log10(fMs),rMTF,'.-k','linewidth',3)
set(gca,'xtick',log10(fMs),'xticklabel',num2cell(fMs))
xlabel('Modulation Frequency (Hz)');ylabel('Firing Rate (Hz)')
title('Rate Modulation Transfer Function');shg
end
end