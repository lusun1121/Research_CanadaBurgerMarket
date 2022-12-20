% Specification 3: First Stage Estimation

% Settings
specification = 3;                    % Specification
outputdir = 'results';                % Output directory
myseed = 10177101;                    % Random number generator seed.
tol = 1e-5;                           % Optimization tolerance.
maxiter = 50000;                      % Maximum number of iterations.
maxeval = 2*maxiter;                  % Maximum number of functional evaluations.
numtheta = 91;
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
% Search near previously estimated parameters for a different specification.
numparticles = 1000;
numstart = 24;
tol = 1e-3;                         % Optimization tolerance.
maxiter = 5000;                     % Maximum number of iterations.
theta_spec12 = [ -0.152486969082878; 0.075075306143974; -0.032963069228375; 0.277389977407132; -0.143824484768414; -0.154579503168450; 0.049657083924141; 0.071839816612254; -0.059796566317599; -0.034050021383283; -0.015329176153995; -0.133410614417299; 0.044937317391017; 0.161510567098635; 0.010427484303121; -0.016703810965091; 0.156063434942925; -0.003422026392538; -0.000902295374939; 0.005076886966493; 0; 0; 0; 0; 0; -5.082670036419369; -4.280381717864108; -3.140533661614120; 1.123034276625564; 2.136400341892982; 2.764217246526321; 3.848953116445719; 0.045922118990243; 0.055998429796226; 0.811289477287182; 0.179897251818038; 0.124441794073815; -0.022914012744017; 0.007140873351643; 0.202918158404322; 0.084035352612383; 0.000207363552419; 0.000158398156369; -0.000000235664521; 1.178336599069278; 0.000537288505304; 0.000028725510431; -0.000000039147560; 0.848416041886465; 0.000365787753789; 0.000098587376475; -0.000000004598740; 0.964045119433548; 0.000292391179604; 0.000019568610453; -0.000000184293062; 0.716383142061953; -0.000432550768017; 0.000259016780781; -0.000000045273453; 0.599144201518315; 0.110028717955094; 0.414568166755614; -0.118082355033924; 0.937708258763673; 0.720296307372555; 0.025914184112982; 0.173039066566827; 0.046117478972812; -0.014239499878435; 0.152053074236773; 0.154751013428532; 0.429227336584324; -0.521517930277520; -0.034176651065707; 0.391713375060937; -0.048075496771015; 0.652234728456667; 0.002013008209404; 0.001564227133266; -0.175788774033504; -0.310629928249856; -0.238769629064274; -0.112858278182169; -0.184589868728525; -0.582777614155305; 0.053200354020951; 0.239002103313829; 0.089337635646161; 0.034693803216762; 0.463413004971719 ];
o = ones(numstart,1);
idx_rz = [ 21 ];
idx_int = [ 22 23 24 25 ];
theta_start = kron(o, theta_spec12');
theta_start(2:end,idx_rz) = -abs(normrnd(0, 0.5, numstart-1, 1)); % Rival Z
theta_start(2:end,idx_int) = normrnd(0, 0.01, numstart-1, 4); % Z interactions

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
        disp(sprintf('%4d  %s', j, sprintf('%9.4f', fval, theta')));
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
