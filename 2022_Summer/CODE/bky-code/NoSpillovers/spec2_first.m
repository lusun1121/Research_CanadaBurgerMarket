% Specification 2: First Stage Estimation

% Settings
format long;
specification = 2;                    % Specification
outputdir = 'results';                % Output directory
myseed = 10177101;                    % Random number generator seed.
tol = 1e-5;                           % Optimization tolerance.
maxiter = 50000;                      % Maximum number of iterations.
maxeval = 2*maxiter;                  % Maximum number of functional evaluations.
numtheta = 72;
numplayers = 5;

% Load the fast food data in MATLAB matrix format.
% Data is city-time panel, where time is indexed by year.
load canadafastfood

% Scaling of the variables.
stdX = sqrt(var(data));

% Optimization settings
options = optimset('fminsearch');
options = optimset(options, 'TolFun', tol);
options = optimset(options, 'MaxFunEvals', maxeval);
options = optimset(options, 'MaxIter', maxiter);
options = optimset(options, 'Display', 'off');

% First-stage starting values
% Search near previously estimated parameters for a nested specification.
numparticles = 1000;
numstart = 24;
theta_spec11_2 = [ -0.14578715027481021; 0.20993124646748884; -0.14474483943408561; 0.14346668119384437; -0.07879363975504577; -0.13044779100448972; 0.06889628474735854; 0.04443159189162982; -0.06682008236261078; -0.05350736435308902; -0.02382506885710049; -0.13753592678436294; 0.03442066757139270; 0.17151490941156333; 0.01374810026774034; -0.01390358943236060; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; -3.89755907751130337; -3.16555845824473892; -2.40109290040271439; 0.85439589938674065; 1.61664222343809705; 2.08934773734107893; 2.85361276546019038; 0.01290885549137999; 0.00841562121418053; 0.43332185678718793; -0.00316753404371494; 0.22004024733419258; 0.07280666844146377; 0.09575010920050564; 0.37291390455699081; 0.22968650588679468; 0.38183476133496197; 0.12752230628398376; 0.40248430221049547; -0.06124983856933999; 1.17934987870176533; 0.81782070812477858; 0.08081731411467188; 0.21258904671216844; 0.25756250586384921; -0.00756446845702947; 0.13885034530131046; 0.33503181304750385; 0.42526446771001780; -0.35619457640318997; 0.71827671114388592; 0.23025473970870952; 0.41633194429605402; 0.82797157472427907; -0.09134498576128641; 0.00225674916808951; -0.13633912462169656; -0.32886522588054357; -0.19036889164125004; -0.08314432382792214; -0.14683794953680979; -0.44904679129455261; 0.95735132348788365; 0.79011149308749551; 0.08745878792524028; 0.02899970468344129; 0.28213093210131057; ];
o = ones(numstart,1);
idx_int = [ 17 18 19 20 22 23 24 25 ];
idx_rz = [ 21 ];
theta_start = kron(o, theta_spec11_2');
theta_start(2:end,idx_int) = normrnd(0, 0.01, numstart-1, 8); % Z interactions
theta_start(2:end,idx_rz) = -abs(normrnd(0, 0.5, numstart-1, 1)); % Rival Z

% Post-process starting values
theta_start(:,26:32) = sort(theta_start(:,26:32), 2); % Sort cutoffs

% Output file names
basename = sprintf('%s/spec%d_first', outputdir, specification);
logname = sprintf('%s.log', basename);
savename = sprintf('%s.mat', basename);

% Start worker pool
ownpool = 0;
if numstart > 1
    if matlabpool('size') == 0
        matlabpool();
        ownpool = 1;
    end
else
    options = optimset(options, 'Display', 'iter-detailed');
end

% Open log file
diary(logname);

% Reporting
disp('Spillovers First-Stage Estimation');
disp('=================================');
disp('');
disp(sprintf('Number of particles: %d', numparticles));
disp(sprintf('Maximum number of iterations: %d', maxiter));
disp(sprintf('Tolerance: %g', tol));
disp(sprintf('RNG Seed: %d', myseed));
disp(sprintf('Specification: %d', specification));
disp('');

% Allocate storage for results
theta_opt = zeros(numstart, numtheta);
fval_opt = zeros(numstart, 1);
exitflag = zeros(numstart, 1);

% Optimize the objective function using each starting value
if numstart > 1
    disp(sprintf('%4s  %9s  %9s', '#   ', 'Function ', 'Param.   '));
end
parfor j = 1:numstart
    [theta, fval, flag] = fminsearch('loglik', theta_start(j,:)', options, data, stdX, specification, numparticles, myseed);
    if numstart > 1
        disp(sprintf('%4d  %s', j, sprintf('%12.5e', fval, theta')));
    end
    theta_opt(j,:) = theta';
    fval_opt(j) = fval;
    exitflag(j) = flag;
end

% Determine the best starting value
fval_best = min(fval_opt);
idx_best = find(fval_opt == fval_best);
theta_best = theta_opt(idx_best, :)';

% Report optimum
disp('Best functional value:');
fval_best
disp('Best parameter values:');
theta_best

% Save results
save(savename, 'specification', 'numparticles', 'theta_start', 'theta_best', 'fval_best', 'idx_best', 'fval_opt', 'theta_opt', 'exitflag');

% Close log file
diary off

% Close worker pool
if numstart > 1 && ownpool > 0
    matlabpool('close');
end

exit;
