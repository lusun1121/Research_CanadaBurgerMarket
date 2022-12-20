% Dynamic Spillovers in the Retail Industry.       %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang. %
% Main program for bootstrap standard errors.      %
% August 31, 2013.                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First and second stage settings
outputdir = 'results';                % Output directory
maxiter = 5000;                       % Maximum number of iterations.
maxeval = 2*maxiter;                  % Maximum number of functional evaluations.
numpaths = 250;                       % Number of paths to simulate.
discount = 0.95;                      % Discount rate.
numperturb = 3000;                    % Number of perturbations.
numperiods = 100;                     % Number of paths to simulate.
tol = 1e-3;                           % Optimization tolerance.
sim_seed = 949751012;                 % Set the seed value for simulations.
ll_seed = 10177101;                   % Set the seed for particle filtering.

% Perturbation standard deviations
sigma = zeros(34, numperturb);

sigma(1:27,     1:250)  = 0.0001 * ones(27, 250);
sigma(1:27,   251:500)  = 0.0010 * ones(27, 250);
sigma(1:27,   501:750)  = 0.0100 * ones(27, 250);
sigma(1:27,   751:1000) = 0.1000 * ones(27, 250);
sigma(28:34,    1:1000) = 0.0100 * ones(7,  1000);

sigma(1:27,  1001:1250) = 0.0001 * ones(27, 250);
sigma(1:27,  1251:1500) = 0.0010 * ones(27, 250);
sigma(1:27,  1501:1750) = 0.0100 * ones(27, 250);
sigma(1:27,  1751:2000) = 0.1000 * ones(27, 250);
sigma(28:34, 1001:2000) = 0.1000 * ones(7,  1000);

sigma(1:27,  2001:2250) = 0.0001 * ones(27, 250);
sigma(1:27,  2251:2500) = 0.0010 * ones(27, 250);
sigma(1:27,  2501:2750) = 0.0100 * ones(27, 250);
sigma(1:27,  2751:3000) = 0.1000 * ones(27, 250);
sigma(28:34, 2001:3000) = 1.0000 * ones(7,  1000);

% Load the fast food data in MATLAB matrix format.
% Data is city-time panel, where time is indexed by year.
load canadafastfood

% Scaling of the variables.
stdX = sqrt(var(data));

% First-stage optimization settings
options = optimset('fminsearch');
options = optimset(options,'LargeScale','Off');
options = optimset(options,'TolFun',tol);
options = optimset(options,'Display','iter-detailed');
options = optimset(options,'MaxFunEvals',maxeval);
options = optimset(options,'HessUpdate','bfgs');
options = optimset(options,'FunValCheck','on');
options = optimset(options,'InitialHessType','scaled-identity');
options = optimset(options,'MaxIter',maxiter);
options = optimset(options,'Diagnostics','on');

% Nonlinear least squares options
lb = [];
ub = [];
%optionsS = optimoptions('lsqnonlin', 'Jacobian', 'off', 'Display', 'iter-detailed');
optionsS = optimset('lsqnonlin');
optionsS = optimset('Jacobian', 'off');
optionsS = optimset('Display', 'iter-detailed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bootstrap settings.
if ~exist('numbootstrap')
    numbootstrap = 100;                 % Number of bootstrap samples.
end
numsamples = 31;
blocklength = 36;

% Determine which specifications to run in first stage.
if ~exist('specification')
    specification = 9;
end
if ~exist('numparticles')
    numparticles = 1000;
end
if ~exist('numthreads')
    numthreads = 1;
end

% Set 'loadfile' to be the filename containing theta_best and thetaS_best
if exist('loadfile')
    load(loadfile)
else
    error('Variable `loadfile` must be set before running main_bootstrap.')
end

% Set base seed for data subsampling
if ~exist('baseseed')
    baseseed = 645456;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output file names
basename = sprintf('%s/spec%d_bootstrap', outputdir, specification);
logname = sprintf('%s.log', basename);
savename = sprintf('%s.mat', basename);

% Number of parameters in first and second stages.
numpar = size(theta_best, 1);
thetaS_best = thetaS_best(:,1);
numparS = size(thetaS_best, 1);

% Create matrices that store parameter estimates with each bootstrap iteration.
theta_boot = zeros(numpar,numbootstrap);
thetaS_boot = zeros(numparS,numbootstrap);

% Open log file
diary(logname);

if numthreads > 1
    matlabpool(numthreads);
    options = optimset(options,'Display','off');
    %optionsS = optimoptions(optionsS, 'Display', 'off');
    optionsS = optimset(optionsS, 'Display', 'off');
end

parfor i = 1:numbootstrap

    % Report progress
    disp(sprintf('Bootstrap replication %d...', i));

    % Subsample the data
    years = 36;
    nummarkets = length(data) / years;
    myseed = baseseed + 27437*i;                            % Set the seed value for randomly drawn subsamples.
    rng(myseed);                                            % Set the seed.
    dataS = datasubsampblock(numsamples,blocklength,data);                   % Parse out data associated with sampled indices.

    % Obtain estimates for SUR regression with subsample.
    [rhoX, sigmaX] = sur(dataS);

    % First stage estimation using subsample.
    [theta, fval, flag] = fminsearch('loglik', theta_best, options, dataS, stdX, specification, numparticles, ll_seed);
    theta_boot(:,i) = theta;

    % Report progress
    disp(sprintf('Replication %4d, first stage:', i));
    disp(sprintf('Log-likelihood: %9.4f', fval));
    disp('Estimates:');
    disp(theta');

    % Second stage estimation using subsample.
    [APIstar, APIperturb] = genperturb(dataS, stdX, theta, rhoX, sigmaX, specification, numpaths, numperiods, discount, numperturb, sigma, sim_seed);
    [thetaS, fvalS, resid, flag] = lsqnonlin('mindist_nlls', thetaS_best, lb, ub, optionsS, APIstar, APIperturb, specification);
    thetaS_boot(:,i) = thetaS;

    % Report progress
    disp(sprintf('Replication %4d, second stage:', i));
    disp(sprintf('Minimum distance function: %9.4f', fvalS));
    disp('Estimates:');
    disp(thetaS');

end

% Close pool
if numthreads > 1
    matlabpool('close');
end

% Store the first and second stage standard errors from bootstrap.
theta_se = sqrt(var(theta_boot'));
thetaS_se = sqrt(var(thetaS_boot'));

% Report optimum
disp('theta_se:');
theta_se
disp('thetaS_se:');
thetaS_se

% Save results
save(savename, 'specification', 'theta_best', 'thetaS_best', 'theta_boot', 'thetaS_boot', 'theta_se', 'thetaS_se');

% Close diary
diary off
