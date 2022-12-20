% Specification 1: First Stage Estimation

% Settings
specification = 1;                    % Specification
%outputdir = 'results';                % Output directory
myseed = 10177101;                    % Random number generator seed.
tol = 1e-5;                           % Optimization tolerance.
maxiter = 50000;                      % Maximum number of iterations.
maxeval = 2*maxiter;                  % Maximum number of functional evaluations.
numtheta = 57;
numparticles = 0;

% Load the fast food data in MATLAB matrix format.
% Data is city-time panel, where time is indexed by year.
load canadafastfood

% Scaling of the variables.
stdX = sqrt(var(data));

% Optimization settings
options = optimoptions('fminunc');
options = optimoptions(options, 'Algorithm', 'quasi-newton');
options = optimoptions(options, 'MaxFunEvals', maxeval);
options = optimoptions(options, 'MaxIter', maxiter);
options = optimoptions(options, 'TolFun', tol);
options = optimoptions(options, 'Display', 'off');
%options = optimoptions(options, 'Display', 'iter-detailed');

% Start worker pool
%threads = matlabpool('size');
delete(gcp('nocreate'));
%threads = parpool('threads');
% if threads < 1
%     matlabpool('open');
%     threads = matlabpool('size');
% end

% Start with ordered probit estimates from Stata
theta_stata = [ -.1415898; .2081491; -.1398425; .1479199; -.1394528; -.1251621; .055328; .0486884; -.064598; -.0476361; -.0232923; -.1323229; .0354094; .1709373; .0139131; -.0146488; -3.466304; -2.84666; -2.181145; .8185702; 1.52743; 1.948664; 2.637176 ; -.0944945; -.0658257; .6408552; -.0045667; .1063918; .3671427; -.0622312; 1.354735; .920912; .0749279; .2101155; .2314288; -.0088738; .1782534; .3253802; .4015002; -.3360856; .7431793; .2260503; .4473324; .9488305; -.0790994; .0021278; -.1506345; -.2964931; -.1582441; -.0824332; -.1684278; -.445869; .9219383; .8253207; .0713067; .0269466; .3192049 ];

% Starting values and settings
numstart = 1;
theta_start = ones(numstart, numtheta);
theta_start(1,:) = theta_stata';
options = optimoptions(options, 'Display', 'iter-detailed');

% Post-process starting values
my = sort(theta_start(:,17:23), 2); % Sort cutoffs

% Output file names
%basename = sprintf('%s/spec%d_first', outputdir, specification);
basename = sprintf('spec%d_first',specification);
%logname = sprintf('%s.log', basename);
savename = sprintf('%s.mat', basename);

% Open log file
%diary(logname);

% Reporting
% disp('Spillovers First-Stage Estimation');
% disp('=================================');
% disp('');
% disp(sprintf('Number of particles: %d', numparticles));
% disp(sprintf('Maximum number of iterations: %d', maxiter));
% disp(sprintf('Tolerance: %g', tol));
% disp(sprintf('RNG Seed: %d', myseed));
% disp(sprintf('Specification: %d', specification));
% disp('');
sprintf('Spillovers First-Stage Estimation');
sprintf('=================================');
sprintf('');
sprintf('Number of particles: %d', numparticles);
sprintf('Maximum number of iterations: %d', maxiter);
sprintf('Tolerance: %g', tol);
sprintf('RNG Seed: %d', myseed);
sprintf('Specification: %d', specification);
sprintf('');

% Allocate storage for results
theta_opt = zeros(numstart, numtheta);
fval_opt = zeros(numstart, 1);
exitflag = zeros(numstart, 1);

% Optimize the objective function using each starting value
% disp(sprintf('%4s  %9s  %9s', '#   ', 'Function ', 'Param.   '));
sprintf('%4s  %9s  %9s', '#   ', 'Function ', 'Param.   ')
for j = 1:numstart
%parfor j = 1:numstart
    [theta, fval, flag] = fminunc('loglik', theta_start(j,:)', options, data, stdX, specification, numparticles, myseed);
%     disp(sprintf('%4d  %s', j, sprintf('%9.4f', fval, theta')));
    sprintf('%4d  %s', j, sprintf('%9.4f', fval, theta'))

    theta_opt(j,:) = theta';
    fval_opt(j) = fval;
    exitflag(j) = flag;
end

% Determine the best starting value
fval_best = min(fval_opt);
idx_best = find(fval_opt == fval_best);
theta_best = theta_opt(idx_best, :)';

% % Report optimum
% disp('Best functional value:');
% fval_best
% disp('Best parameter values:');
% theta_best
sprintf('Best functional value:');
fval_best
sprintf('Best parameter values:');
theta_best

% Save results
save(savename, 'specification', 'numparticles', 'theta_start', 'theta_best', 'fval_best', 'idx_best', 'fval_opt', 'theta_opt', 'exitflag');

% Close log file
% diary off

% Close worker pool
%matlabpool('close');
%delete(gcp('nocreate'));
%exit;
