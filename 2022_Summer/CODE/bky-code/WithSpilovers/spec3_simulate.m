% Specification 3: Forward Simulation

% Settings
specification = 3;                    % Specification
discount = 0.95;                      % Discount rate.
myseed = 949751012;                   % Set the seed value for simulations.

% Simulation settings
numpaths = 250;                   % Number of paths to simulate.
numperiods = 100;                 % Number of paths to simulate.
numperturb = 3000;                % Number of perturbations.
outputdir = 'results';            % Output directory

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

% Load first stage results
load(sprintf('%s/spec%d_first.mat', outputdir, specification));

% Output file names
savename = sprintf('%s/spec%d_simulate.mat', outputdir, specification);

% Load dataset.
load canadafastfood

% Scaling of the variables.
stdX = sqrt(var(data));

% Obtain SUR estimation results.
[rhoX, sigmaX] = sur(data);

% Start worker pool
threads = matlabpool('size');
if threads < 1
    matlabpool('open');
    threads = matlabpool('size');
end

% Precalculate simulations using first-stage estimates
[APIstar, APIperturb] = genperturb(data, stdX, theta_best, rhoX, sigmaX, specification, numpaths, numperiods, discount, numperturb, sigma, myseed);

% Save results
save(savename, 'specification', 'theta_best', 'rhoX', 'sigmaX', 'APIstar', 'APIperturb');

% Close worker pool
matlabpool('close');

exit;
