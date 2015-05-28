
% rateFromSorted: Counts # of spikes for data stored in cell-upon-cell
% format. Makes it easier to remove neurons from analysis based on sparseness
% - MJRunfeldt April 2015

% INPUT = spiketimes {neuron} {stimulus} {trial}. Temporal resolution
% irrelevant because output is a count, not a rate.
% OUTPUT: (1) total number of spikes for each neuron
% (2) total number of spikes for each modulation frequency
% (3) Average number of spikes per neuron 

function [spNrn,spFM,nrnMean] = rateFromSorted(input)
% % input= spikes.JPsort ; % for devel

nNrns = length(input);

spFM = cell(1,length(nNrns));
for n = 1:nNrns
    iNrn = input{n}; % {fm} {trial}
    spTrial = cell(1,length(iNrn));
    for fM = 1:length(iNrn)
        spTrial{fM} = cell2mat(cellfun(@length,iNrn{fM},'un',0)) ; % # of spikes per trial
    end
    spFM{n} = cell2mat(cellfun(@sum,spTrial,'un',0)) ; % # of spikes per fM
end

spNrn = cell2mat(cellfun(@sum,spFM,'un',0)) ; % total # of spikes per neuron
nrnMean = cell2mat(cellfun(@mean,spFM,'un',0)) ; % grand average # of spikes per

end
