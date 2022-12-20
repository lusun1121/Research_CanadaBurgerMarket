% Dynamic Spillovers in the Retail Industry.        %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.  %
% Forward simulation using estimated policies.      %
% October  8, 2014                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function APIstar = genstar(data, stdX, param, rhoX, corrX, specification, numpaths, numperiods, discount, numperturb, myseed)

% Parameters for main settings.
years = 36;
numplayers = 5;                     % Number of players.
nummarkets = length(data) / years;  % Number of markets.
rng(myseed);                        % Set the seed.
ng = 7;                             % Number of cutoffs.
%numtheta = 17 + ng;                 % Number of policy parameters per player.
numtheta = 27 + ng;                 % Number of policy parameters per player.

% Set perturbations to zero
noperturb = zeros(numtheta*numplayers, 1);

% Simulated values using estimated first stage parameters.
disp('Simulating equilibrium policies...');
[APIAWs, APIBKs, APIHARVs, APIMCDs, APIWENDs] = forwardsim(param, noperturb, rhoX, corrX, data, stdX, specification, numpaths, numperiods, discount);

% Stack equilibrium simulation for comparison with alternatives
o = ones(numperturb,1);
APIAWstar = kron(o, APIAWs);
APIBKstar = kron(o, APIBKs);
APIHARVstar = kron(o, APIHARVs);
APIMCDstar = kron(o, APIMCDs);
APIWENDstar = kron(o, APIWENDs);

AO = zeros(size(APIAWstar,1),size(APIAWstar,2));
APIstar = [APIAWstar AO AO AO AO; ...
           AO APIBKstar AO AO AO; ...
           AO AO APIHARVstar AO AO; ...
           AO AO AO APIMCDstar AO; ...
           AO AO AO AO APIWENDstar];
end
