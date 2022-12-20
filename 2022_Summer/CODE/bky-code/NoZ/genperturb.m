% Dynamic Spillovers in the Retail Industry.        %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.  %
% Forward simulation using perturbed policies.      %
% October  8, 2014                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [APIstar, APIperturb] = genperturb(data, stdX, param, rhoX, corrX, specification, numpaths, numperiods, discount, numperturb, sigma, myseed)

% Parameters for main settings.
years = 36;
numplayers = 5;                     % Number of players.
nummarkets = length(data) / years;  % Number of markets.
rng(myseed);                        % Set the seed.
ng = 7;                             % Number of cutoffs.
numtheta = 27 + ng;                 % Number of policy parameters per player.

% Firm-specific perturbations of ordered probit parameters
mu = zeros(numtheta, numperturb);
perturb = normrnd(mu, sigma); % Single-player perturbations.
PO = zeros(numtheta, numperturb);
perturbAW = [ perturb; PO; PO; PO; PO; ];
perturbBK = [ PO; perturb; PO; PO; PO; ];
perturbHARV = [ PO; PO; perturb; PO; PO; ];
perturbMCD = [ PO; PO; PO; perturb; PO; ];
perturbWEND = [ PO; PO; PO; PO; perturb; ];

% Initialize storage for simulations
APIAWperturb = [];
APIBKperturb = [];
APIHARVperturb = [];
APIMCDperturb = [];
APIWENDperturb = [];

% Simulated values using estimated first stage parameters.
APIstar = genstar(data, stdX, param, rhoX, corrX, specification, numpaths, numperiods, discount, numperturb, myseed);

parfor i = 1:numperturb
    disp(sprintf('Simulating with alternative policy %d...', i));

    % Simulated values using perturbed first stage parameters for each firm.
    [APIAWpb, APIBK, APIHARV, APIMCD, APIWEND] = forwardsim(param, perturbAW(:,i), rhoX, corrX, data, stdX, specification, numpaths, numperiods, discount);
    [APIAW, APIBKpb, APIHARV, APIMCD, APIWEND] = forwardsim(param, perturbBK(:,i), rhoX, corrX, data, stdX, specification, numpaths, numperiods, discount);
    [APIAW, APIBK, APIHARVpb, APIMCD, APIWEND] = forwardsim(param, perturbHARV(:,i), rhoX, corrX, data, stdX, specification, numpaths, numperiods, discount);
    [APIAW, APIBK, APIHARV, APIMCDpb, APIWEND] = forwardsim(param, perturbMCD(:,i), rhoX, corrX, data, stdX, specification, numpaths, numperiods, discount);
    [APIAW, APIBK, APIHARV, APIMCD, APIWENDpb] = forwardsim(param, perturbWEND(:,i), rhoX, corrX, data, stdX, specification, numpaths, numperiods, discount);

    APIAWperturb = [APIAWperturb; APIAWpb];
    APIBKperturb = [APIBKperturb; APIBKpb];
    APIHARVperturb = [APIHARVperturb; APIHARVpb];
    APIMCDperturb = [APIMCDperturb; APIMCDpb];
    APIWENDperturb = [APIWENDperturb; APIWENDpb];
end

% Store the equilibrium and perturbed discounted sums
AO = zeros(size(APIAWperturb,1),size(APIAWperturb,2));
APIperturb = [APIAWperturb AO AO AO AO; ...
              AO APIBKperturb AO AO AO; ...
              AO AO APIHARVperturb AO AO; ...
              AO AO AO APIMCDperturb AO; ...
              AO AO AO AO APIWENDperturb];
end
