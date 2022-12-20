% Dynamic Spillovers in the Retail Industry.             %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.       %
% Main program for simulation analysis.                  %
% January 6, 2015                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
numpaths = 10000;                   % Number of paths to simulate.
discount = 0.95;                    % Discount rate.
numperiods = 36;                    % Number of years to simulate.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changes to initial conditions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seed the random number generator
rng(myseed);

% Changing the initial conditions for all retailers.
alpha = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
[meanshareinit1] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 0, alpha(1));
[meanshareinit2] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 0, alpha(2));
[meanshareinit3] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 0, alpha(3));
[meanshareinit4] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 0, alpha(4));
[meanshareinit5] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 0, alpha(5));
[meanshareinit6] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 0, alpha(6));
[meanshareinit7] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 0, alpha(7));
[meanshareinit8] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 0, alpha(8));
[meanshareinit9] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 0, alpha(9));
[meanshareinit10] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 0, alpha(10));
[meanshareinit11] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 0, alpha(11));

% Creates table of counterfactual results in LaTeX.
columnLabels = {'0', '1', '2', '3', '4', '5' , '6', '7', '8', '9', '10'};
rowLabels = {'A\&W', 'Burger King', 'Harvey''s', 'McDonald''s', 'Wendy''s', 'HHI'};
results = [meanshareinit1; meanshareinit2; meanshareinit3; meanshareinit4; meanshareinit5; meanshareinit6; meanshareinit7; meanshareinit8; meanshareinit9; meanshareinit10; meanshareinit11];
results = results';
hhi = results(1,:).^2 + results(2,:).^2 + results(3,:).^2 + results(4,:).^2 + results(5,:).^2;
results = [results; hhi];
matrix2latex(results, 'counterfactualresults_initconditions.tex', 'columnLabels', columnLabels, 'rowLabels', rowLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'footnotesize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impact of economic changes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seed the random number generator
rng(myseed);

% Looking at the impact that recession has on McDonald's dominance.
alpha = [-0.1, -0.05, 0.05, 0.1];
[meanshareecon1] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 1, alpha(1));
[meanshareecon2] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 1, alpha(2));
[meanshareecon3] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 1, alpha(3));
[meanshareecon4] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 1, alpha(4));

% Creates table of counterfactual results in LaTeX.
columnLabels = {'-10\%', '-5\%', '5\%', '10\%'};
results = [meanshareecon1; meanshareecon2; meanshareecon3; meanshareecon4];
results = results';
hhi = results(1,:).^2 + results(2,:).^2 + results(3,:).^2 + results(4,:).^2 + results(5,:).^2;
results = [results; hhi];
matrix2latex(results, 'counterfactualresults_econconditions.tex', 'columnLabels', columnLabels, 'rowLabels', rowLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'footnotesize');

% Seed the random number generator
rng(myseed);

% Looking at the impact that minimum wage changes have on McDonald's dominance.
alpha = [-0.1, -0.05, 0.05, 0.1];
[meanshareecon1] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 2, alpha(1));
[meanshareecon2] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 2, alpha(2));
[meanshareecon3] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 2, alpha(3));
[meanshareecon4] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 2, alpha(4));

% Creates table of counterfactual results in LaTeX.
columnLabels = {'-10\%', '-5\%', '5\%', '10\%'};
results = [meanshareecon1; meanshareecon2; meanshareecon3; meanshareecon4];
results = results';
hhi = results(1,:).^2 + results(2,:).^2 + results(3,:).^2 + results(4,:).^2 + results(5,:).^2;
results = [results; hhi];
matrix2latex(results, 'counterfactualresults_wageconditions.tex', 'columnLabels', columnLabels, 'rowLabels', rowLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'footnotesize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endow all firms with Z's from stationary distribution of McDonald's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seed the random number generator
rng(myseed);

alpha = [1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 30.0];
[meansharez1] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 3, alpha(1));
[meansharez2] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 3, alpha(2));
[meansharez3] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 3, alpha(3));
[meansharez4] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 3, alpha(4));
[meansharez5] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 3, alpha(5));
[meansharez6] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 3, alpha(6));
[meansharez7] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 3, alpha(7));
[meansharez8] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 3, alpha(8));

% Creates table of counterfactual results in LaTeX.
columnLabels = {'Baseline', '100\%', '200\%', '300\%', '400\%', '500\%', '1000\%', '2000\%', '3000\%'};
results = [ meanshareinit1; meansharez1; meansharez2; meansharez3; meansharez4; meansharez5; meansharez6; meansharez7; meansharez8; ]; % Use precalculated baseline values (meanshareinit1)
results = results';
hhi = results(1,:).^2 + results(2,:).^2 + results(3,:).^2 + results(4,:).^2 + results(5,:).^2;
results = [results; hhi];
matrix2latex(results, 'counterfactualresults_spillovers.tex', 'columnLabels', columnLabels, 'rowLabels', rowLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'footnotesize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw initial Z's from stationary distributions with rescaled delta's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seed the random number generator
rng(myseed);

alpha = [ 1.0, 0.5, 0.0, -0.5, -1.0 ];
[meanshareforget1] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 4, alpha(1));
[meanshareforget2] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 4, alpha(2));
%[meanshareforget3] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 4, alpha(3));
meanshareforget3 = meanshareinit1; % Use precalculated baseline values
[meanshareforget4] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 4, alpha(4));
[meanshareforget5] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 4, alpha(5));

% Creates table of counterfactual results in LaTeX.
columnLabels = {'$2 \delta_i$', '$1.5 \delta_i$', '$\delta_i$', '$0.5 \delta_i$', '$0$' };
results = [ meanshareforget1; meanshareforget2; meanshareforget3; meanshareforget4; meanshareforget5; ];
results = results';
hhi = results(1,:).^2 + results(2,:).^2 + results(3,:).^2 + results(4,:).^2 + results(5,:).^2;
results = [ results; hhi ];
matrix2latex(results, 'counterfactualresults_forget.tex', 'columnLabels', columnLabels, 'rowLabels', rowLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'footnotesize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw initial Z's for McDonald's from distribution with scaled delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seed the random number generator
rng(myseed);

alpha = [ 1.0, 0.5, 0.0, -0.5, -1.0 ];
[meanshareforgetmcd1] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 5, alpha(1));
[meanshareforgetmcd2] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 5, alpha(2));
%[meanshareforgetmcd3] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 5, alpha(3));
meanshareforgetmcd3 = meanshareinit1; % Use precalculated baseline values
[meanshareforgetmcd4] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 5, alpha(4));
[meanshareforgetmcd5] = forwardsimrobust(theta1_spillovers, rhoX, corrX, data, stdX, spillovers_spec1, numpaths, numperiods, discount, theta2_spillovers, spillovers_spec2, 5, alpha(5));

% Creates table of counterfactual results in LaTeX.
columnLabels = {'$2 \delta_i$', '$1.5 \delta_i$', '$\delta_i$', '$0.5 \delta_i$', '$0$' };
results = [ meanshareforgetmcd1; meanshareforgetmcd2; meanshareforgetmcd3; meanshareforgetmcd4; meanshareforgetmcd5; ];
results = results';
hhi = results(1,:).^2 + results(2,:).^2 + results(3,:).^2 + results(4,:).^2 + results(5,:).^2;
results = [ results; hhi ];
matrix2latex(results, 'counterfactualresults_forgetmcd.tex', 'columnLabels', columnLabels, 'rowLabels', rowLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'footnotesize');

exit;
