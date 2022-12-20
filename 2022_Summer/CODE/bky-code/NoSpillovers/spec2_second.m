% Specification 2: Second Stage Estimation

% Settings
specification = 2;

% Load the fast food data in MATLAB matrix format.
% Data is city-time panel, where time is indexed by year.
load canadafastfood

% Scaling of the variables.
stdX = sqrt(var(data));

% Obtain SUR estimation results.
[rhoX, sigmaX] = sur(data);

% Nonlinear least squares options
lb = [];
ub = [];
options = optimoptions('lsqnonlin', 'Jacobian', 'off', 'Display', 'off');
%options = optimoptions(options, 'Display', 'iter-detailed');
%options = optimoptions(options, 'DerivativeCheck', 'on');

% Second-stage settings
outputdir = 'results';                % Output directory
numthetaS = 12;
numstart = 12;
thetaS_start = normrnd(0, 1, numstart, numthetaS);
thetaS_start(1,:) = zeros(1, numthetaS);

% Output file names
basename = sprintf('%s/spec%d_second', outputdir, specification);
logname = sprintf('%s.log', basename);
savename = sprintf('%s.mat', basename);

% Load simulations
load(sprintf('%s/spec%d_simulate.mat', outputdir, specification));

% Start worker pool
if numstart > 1
    threads = matlabpool('size');
    if threads < 1
        matlabpool('open');
        threads = matlabpool('size');
    end
else
    options = optimoptions(options, 'Display', 'iter-detailed');
end

% Open log file
diary(logname);

% Reporting
disp('Spillovers Second-Stage Estimation');
disp('==================================');
disp('');
disp(sprintf('Specification: %d', specification));
if numstart == 1
    disp('Starting values:');
    disp(thetaS_start);
end
disp('');

% Allocate storage for results
thetaS_opt = zeros(numstart, numthetaS);
fvalS_opt = zeros(numstart, 1);
exitflagS = zeros(numstart, 1);

% Optimize the objective function using each starting value
if numstart > 1
    disp(sprintf('%4s  %9s  %9s', '#   ', 'Function ', 'Param.   '));
end
parfor j = 1:numstart
    [thetaS, fvalS, resid, flag] = lsqnonlin('mindist_nlls', thetaS_start(j,:)', lb, ub, options, APIstar, APIperturb, specification);
    if numstart > 1
        disp(sprintf('%4d  %s', j, sprintf('%9.4f', fvalS, thetaS')));
    end
    thetaS_opt(j,:) = thetaS';
    fvalS_opt(j) = fvalS;
    exitflagS(j) = flag;
end

% Determine the best starting value
fvalS_best = min(fvalS_opt);
idxS_best = find(fvalS_opt==fvalS_best);
thetaS_best = thetaS_opt(idxS_best, :)';

% Report optimum
disp('Best functional value:');
fvalS_best
disp('Best parameter values:');
thetaS_best

% Save results
save(savename, 'specification', 'theta_best', 'thetaS_start', 'thetaS_opt', 'fvalS_opt', 'exitflagS', 'thetaS_best', 'fvalS_best', 'idxS_best', 'rhoX', 'sigmaX', 'APIstar', 'APIperturb');

% Close diary
diary off

% Close worker pool
if numstart > 1
    matlabpool('close');
end

exit;
