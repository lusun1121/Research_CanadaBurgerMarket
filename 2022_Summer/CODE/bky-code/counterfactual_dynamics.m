% Dynamic Spillovers in the Retail Industry.                                      %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.                                %
% Main program for simulation analysis of number of stores and Z.                 %
% January 6, 2015                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the data and estimated parameters from previous steps.
load canadafastfood                 % Data.

% Load the estimated first-stage parameters from different specifications.
load('results/spec1_param.mat')
theta1_BBL = theta_best;            % Baseline BBL first stage.
theta2_BBL = thetaS_best;           % Baseline BBL second stage.
BBL_spec1 = 1;
BBL_spec2 = 1;
load('results/spec3_param.mat')
theta1_spillovers = theta_best;     % BBL with Z (spillovers) first stage.
theta2_spillovers = thetaS_best;    % BBL with Z (spillovers) second stage.
spillovers_spec1 = 3;
spillovers_spec2 = 3;
load('results/spec2_param.mat')
theta1_nospillovers = theta_best;   % BBL with Z (no spillovers) first stage.
theta2_nospillovers = thetaS_best;  % BBL with Z (no spillovers) second stage.
nospillovers_spec1 = 2;
nospillovers_spec2 = 2;

% Main settings
numpaths = 1000;                    % Number of paths to simulate.
discount = 0.95;                    % Discount rate.
numperiods = 36;                    % Number of periods to simulate.
myseed = 18740192;                  % Common seed for reproducibility

% Obtain estimates for SUR regression.
[rhoX, corrX] = sur(data);

% Scaling of the variables.
stdX = sqrt(var(data));

% Parse out relevant data for model fit evaluation.
time = data(:,1);           % Year.
numAW = data(:,2);          % Number of A & W outlets in city.
numBK = data(:,3);          % Number of Burger King outlets in city.
numHARV = data(:,4);        % Number of Harvey's outlets in city.
numMCD = data(:,5);         % Number of McDonald's outlets in city.
numWEND = data(:,6);        % Number of Wendy's outlets in city.
cityindex = data(:,20);     % Index for city.

% Obtain the simulated N's, PI's, and Z's for full sample.
rng(myseed);
[Nstore_BBL, Nbias_BBL, Nmse_BBL, Zstore_BBL, PIstore_BBL] = forwardsiminfo(theta1_BBL, rhoX, corrX, data, stdX, BBL_spec1, numpaths, numperiods, discount, theta2_BBL, BBL_spec2, 0);

rng(myseed);
[Nstore_spillovers, Nbias_spillovers, Nmse_spillovers, Zstore_spillovers, PIstore_spillovers] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 0);

rng(myseed);
[Nstore_nospillovers, Nbias_nospillovers, Nmse_nospillovers, Zstore_spillovers, PIstore_nospillovers] = forwardsiminfo(theta1_nospillovers, rhoX, corrX, data, stdX, nospillovers_spec1, numpaths, numperiods, discount, theta2_nospillovers, nospillovers_spec2, 0);

model_fit_data = [cityindex, time, numAW, numBK, numHARV, numMCD, numWEND, Nstore_BBL, Nstore_spillovers, Nstore_nospillovers, Nbias_BBL, Nbias_spillovers, Nbias_nospillovers, Nmse_BBL, Nmse_spillovers, Nmse_nospillovers];

csvwrite('model_fit_data.csv',model_fit_data);

% Obtain the simulated PI's.
market_share_data = [cityindex, time, numAW, numBK, numHARV, numMCD, numWEND, PIstore_BBL, PIstore_spillovers, PIstore_nospillovers];

csvwrite('market_share_data.csv', market_share_data);

% Save the simulated Z's.
market_share_Z_data = [cityindex, time, numAW, numBK, numHARV, numMCD, numWEND, Zstore_spillovers];
csvwrite('market_share_Z_data.csv',market_share_Z_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impulse response analysis: one unit increase in Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Impulse response analysis for A&W, increase in Z
rng(myseed);
[Nstore_IRI_AW, Nbias_IRI_AW, Nmse_IRI_AW, Zstore_IRI_AW, PIstore_IRI_AW] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 1);

% Impulse response analysis for Burger King, increase in Z
rng(myseed);
[Nstore_IRI_BK, Nbias_IRI_BK, Nmse_IRI_BK, Zstore_IRI_BK, PIstore_IRI_BK] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 2);

% Impulse response analysis for Harvey's, increase in Z
rng(myseed);
[Nstore_IRI_HARV, Nbias_IRI_HARV, Nmse_IRI_HARV, Zstore_IRI_HARV, PIstore_IRI_HARV] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 3);

% Impulse response analysis for McDonald's, increase in Z
rng(myseed);
[Nstore_IRI_MCD, Nbias_IRI_MCD, Nmse_IRI_MCD, Zstore_IRI_MCD, PIstore_IRI_MCD] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 4);

% Impulse response analysis for Wendy's, increase in Z
rng(myseed);
[Nstore_IRI_WEND, Nbias_IRI_WEND, Nmse_IRI_WEND, Zstore_IRI_WEND, PIstore_IRI_WEND] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 5);

% Store originally simulated data and impulse response simulations
impulse_response_Z_data = [ cityindex, time, numAW, numBK, numHARV, numMCD, numWEND, Nstore_spillovers, Nstore_IRI_AW, Nstore_IRI_BK, Nstore_IRI_HARV, Nstore_IRI_MCD, Nstore_IRI_WEND, PIstore_spillovers, PIstore_IRI_AW, PIstore_IRI_BK, PIstore_IRI_HARV, PIstore_IRI_MCD, PIstore_IRI_WEND ];
csvwrite('impulse_response_Z_data.csv',impulse_response_Z_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impulse response analysis: one unit decrease in Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Impulse response analysis for A&W, decrease in Z
rng(myseed);
[Nstore_IRI_AW, Nbias_IRI_AW, Nmse_IRI_AW, Zstore_IRI_AW, PIstore_IRI_AW] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, -1);

% Impulse response analysis for Burger King, decrease in Z
rng(myseed);
[Nstore_IRI_BK, Nbias_IRI_BK, Nmse_IRI_BK, Zstore_IRI_BK, PIstore_IRI_BK] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, -2);

% Impulse response analysis for Harvey's, decrease in Z
rng(myseed);
[Nstore_IRI_HARV, Nbias_IRI_HARV, Nmse_IRI_HARV, Zstore_IRI_HARV, PIstore_IRI_HARV] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, -3);

% Impulse response analysis for McDonald's, decrease in Z
rng(myseed);
[Nstore_IRI_MCD, Nbias_IRI_MCD, Nmse_IRI_MCD, Zstore_IRI_MCD, PIstore_IRI_MCD] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, -4);

% Impulse response analysis for Wendy's, decrease in Z
rng(myseed);
[Nstore_IRI_WEND, Nbias_IRI_WEND, Nmse_IRI_WEND, Zstore_IRI_WEND, PIstore_IRI_WEND] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, -5);

% Store originally simulated data and impulse response simulations
impulse_response_ZN_data = [ cityindex, time, numAW, numBK, numHARV, numMCD, numWEND, Nstore_spillovers, Nstore_IRI_AW, Nstore_IRI_BK, Nstore_IRI_HARV, Nstore_IRI_MCD, Nstore_IRI_WEND, PIstore_spillovers, PIstore_IRI_AW, PIstore_IRI_BK, PIstore_IRI_HARV, PIstore_IRI_MCD, PIstore_IRI_WEND ];
csvwrite('impulse_response_ZN_data.csv',impulse_response_ZN_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impulse response analysis: 20% decrease (D) in number of stores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Impulse response analysis for A&W
rng(myseed);
[Nstore_IRD_AW, Nbias_IRD_AW, Nmse_IRD_AW, Zstore_IRD_AW, PIstore_IRD_AW] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 11);

% Impulse response analysis for Burger King
rng(myseed);
[Nstore_IRD_BK, Nbias_IRD_BK, Nmse_IRD_BK, Zstore_IRD_BK, PIstore_IRD_BK] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 12);

% Impulse response analysis for Harvey's
rng(myseed);
[Nstore_IRD_HARV, Nbias_IRD_HARV, Nmse_IRD_HARV, Zstore_IRD_HARV, PIstore_IRD_HARV] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 13);

% Impulse response analysis for McDonald's
rng(myseed);
[Nstore_IRD_MCD, Nbias_IRD_MCD, Nmse_IRD_MCD, Zstore_IRD_MCD, PIstore_IRD_MCD] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 14);

% Impulse response analysis for Wendy's
rng(myseed);
[Nstore_IRD_WEND, Nbias_IRD_WEND, Nmse_IRD_WEND, Zstore_IRD_WEND, PIstore_IRD_WEND] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 15);

% Store originally simulated data and impulse response simulations
impulse_response_D_data = [ cityindex, time, numAW, numBK, numHARV, numMCD, numWEND, Nstore_spillovers, Nstore_IRD_AW, Nstore_IRD_BK, Nstore_IRD_HARV, Nstore_IRD_MCD, Nstore_IRD_WEND, PIstore_spillovers, PIstore_IRD_AW, PIstore_IRD_BK, PIstore_IRD_HARV, PIstore_IRD_MCD, PIstore_IRD_WEND ];
csvwrite('impulse_response_D_data.csv',impulse_response_D_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impulse response analysis: 20% increase (I) in number of stores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Impulse response analysis for A&W
rng(myseed);
[Nstore_IRI_AW, Nbias_IRI_AW, Nmse_IRI_AW, Zstore_IRI_AW, PIstore_IRI_AW] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 21);

% Impulse response analysis for Burger King
rng(myseed);
[Nstore_IRI_BK, Nbias_IRI_BK, Nmse_IRI_BK, Zstore_IRI_BK, PIstore_IRI_BK] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 22);

% Impulse response analysis for Harvey's
rng(myseed);
[Nstore_IRI_HARV, Nbias_IRI_HARV, Nmse_IRI_HARV, Zstore_IRI_HARV, PIstore_IRI_HARV] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 23);

% Impulse response analysis for McDonald's
rng(myseed);
[Nstore_IRI_MCD, Nbias_IRI_MCD, Nmse_IRI_MCD, Zstore_IRI_MCD, PIstore_IRI_MCD] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 24);

% Impulse response analysis for Wendy's
rng(myseed);
[Nstore_IRI_WEND, Nbias_IRI_WEND, Nmse_IRI_WEND, Zstore_IRI_WEND, PIstore_IRI_WEND] = forwardsiminfo(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 25);

% Store originally simulated data and impulse response simulations
impulse_response_I_data = [ cityindex, time, numAW, numBK, numHARV, numMCD, numWEND, Nstore_spillovers, Nstore_IRI_AW, Nstore_IRI_BK, Nstore_IRI_HARV, Nstore_IRI_MCD, Nstore_IRI_WEND, PIstore_spillovers, PIstore_IRI_AW, PIstore_IRI_BK, PIstore_IRI_HARV, PIstore_IRI_MCD, PIstore_IRI_WEND ];
csvwrite('impulse_response_I_data.csv',impulse_response_I_data);

exit;
