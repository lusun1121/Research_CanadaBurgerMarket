% Dynamic Spillovers in the Retail Industry.                                      %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.                                %
% Subsample data for block bootstrap.                                             %
% September 9, 2013.                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function datasubsamp = datasubsampblock(numsamples, blocklength, data)

years = 36;
yearset = [1:1:years-blocklength+1];
nummarkets = length(data) / years;
markets = [1:1:nummarkets];
numblocks = years / blocklength;

% Draw market indices
marketsubset = randsample(markets, numsamples, true);

% Draw time indices
timesubset = randsample(yearset, numblocks, true);

% Initialize data matrix
datasubsamp = [];

% For each market, sample time periods
for market = marketsubset
    for t = timesubset
        dataS = data((market-1)*years+t:(market-1)*years+t+blocklength-1,:);
        datasubsamp = [datasubsamp; dataS];
    end
end

end
