% Dynamic Spillovers in the Retail Industry.            %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.      %
% Log-likelihood function.                              %
% July 22, 2013.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sll = loglik(param, data, stdX, specification, numparticles, myseed)

% Settings for simulation.
years = 36;
nummarkets = length(data)/years;
players = 5;
rng(myseed);                        % Set the seed.
index = [1:1:numparticles];         % Index for each particle simulation draw for a given market-time.

% Parsing out the parameters.
[theta, mu, delta, beta, eta] = getparam1(param, specification);

% Initialize matrix which stores log-likelhood for each market-time.
ll = zeros(years,nummarkets);

% Loop over all markets.
for market = 1:nummarkets

    % Parse out the data based on market index.
    dataN = data((market-1)*years+1:market*years,:);

    % Parse out the data.
    time = dataN(:,1);           % Year.
    lagnumAW = dataN(:,7);       % Lagged number of A & W outlets in city.
    lagnumBK = dataN(:,8);       % Lagged number of Burger King outlets in city.
    lagnumHARV = dataN(:,9);     % Lagged number of Harvey's outlets in city.
    lagnumMCD = dataN(:,10);     % Lagged number of McDonald's outlets in city.
    lagnumWEND = dataN(:,11);    % Lagged number of Wendy's outlets in city.
    population = dataN(:,12);    % Population.
    income = dataN(:,13);        % Income.
    value = dataN(:,14);         % Property value.
    aAW = dataN(:,15);           % Change in A & W outlets.
    aBK = dataN(:,16);           % Change in Burger King outlets.
    aHARV = dataN(:,17);         % Change in Harvey's outlets.
    aMCD = dataN(:,18);          % Change in McDonald's outlets.
    aWEND = dataN(:,19);         % Change in Wendy's outlets.
    cityindex = dataN(:,20);     % Index for city.
    timeindex = dataN(:,21);     % Index for time.

    % Exogenous shock to city's exposure.
    greycup = dataN(:,45);              % City hosting CFL Grey Cup.
    regulation = dataN(:,46);           % Introduction of smoking regulation.
    minwage = dataN(:,47);              % Minimum wage.
    lagminwage = dataN(:,48);           % Lagged minimum wage.

    % Lagged actions.
    lagaAW = dataN(:,49);
    lagaBK = dataN(:,50);
    lagaHARV = dataN(:,51);
    lagaMCD = dataN(:,52);
    lagaWEND = dataN(:,53);

    if numparticles > 0
        % Allocate storage for particles for each market
        randparticles = zeros(numparticles*years,players);

        % Store simulation draws used for propogating particles
        randepsilons = normrnd(0,1,numparticles*years,players);
    else
        % Include firm fixed effects via particles
        randparticles = ones(years, 1) * mu' + eta(cityindex) * ones(1, players);
        randepsilons = zeros(years,players);
    end

    % Initialize vector of log-likelihoods which will be stored for future use.
    lmjointP = ones(years,1);

    % Loop over all 36 years.
    for t = 1:years

        if numparticles > 0
            o = ones(numparticles,1);
        else
            o = ones(1,1);
        end

        lagnum = [lagnumAW(t), lagnumBK(t), lagnumHARV(t), lagnumMCD(t), lagnumWEND(t)];
        lagN = kron(lagnum,o);

        if t == 1
            if numparticles > 0
                % Draw Z's from initial distribution in initial time period.
                randparticles(1:numparticles,:) = normrnd(0,1,numparticles,players);
            end
        else
            % Transition Z's in subsequent time periods.
            if numparticles > 0
                lagact = [lagaAW(t), lagaBK(t), lagaHARV(t), lagaMCD(t), lagaWEND(t)];
                lagA = kron(lagact,o);
                Cind = kron(cityindex(t),o);
                re = randepsilons((t-1)*numparticles+1:t*numparticles,:);
                Zhat  = draw_Z(lagZ, lagA, lagN, Cind, re, mu, delta, beta, eta, specification);
                randparticles((t-1)*numparticles+1:t*numparticles,:) = [Zhat(:,1), Zhat(:,2), Zhat(:,3), Zhat(:,4), Zhat(:,5)];
            end
        end

        % Store current swarm
        if numparticles > 0
            rp = randparticles((t-1)*numparticles+1:t*numparticles,:);
        else
            rp = randparticles(t,:);
        end

        % Current actions
        y = [aAW(t); aBK(t); aHARV(t); aMCD(t); aWEND(t)];

        % Rescale variables
        pop = population(t) / stdX(12);
        inc = income(t) / stdX(13);
        val = value(t) / stdX(14);
        host = greycup(t);
        smoke = regulation(t);
        minw = minwage(t) / stdX(47);

        % Stack the state variables for each simulation draw.
        X = [ pop, inc, val, host, smoke, minw, ...
              pop * pop, pop * inc, pop * val, pop * minw, ...
              inc * inc, inc * val, inc * minw, ...
              val * val, val * minw, ...
              minw * minw ];
        X_interact = [ 1, pop, inc, val, minw ];
        Xo = kron(X,o);
        X_AW = [Xo, rp(:,1) * X_interact, mean(rp(:,[2,3,4,5]), 2) * X_interact, lagN(:,1) ];
        X_BK = [Xo, rp(:,2) * X_interact, mean(rp(:,[1,3,4,5]), 2) * X_interact, lagN(:,2) ];
        X_HARV = [Xo, rp(:,3) * X_interact, mean(rp(:,[1,2,4,5]), 2) * X_interact, lagN(:,3) ];
        X_MCD = [Xo, rp(:,4) * X_interact, mean(rp(:,[1,2,3,5]), 2) * X_interact, lagN(:,4) ];
        X_WEND = [Xo, rp(:,5) * X_interact, mean(rp(:,[1,2,3,4]), 2) * X_interact, lagN(:,5) ];

        % Stack the state variables for each simulation draw.
        XM = zeros(size(X_AW,1), size(X_AW,2));
        M = [X_AW XM XM XM XM; XM X_BK XM XM XM; XM XM X_HARV XM XM; XM XM XM X_MCD XM; XM XM XM XM X_WEND];

        % Obtain the probabilities for each simulation draw.
        P = probZ(theta,y,M);
        if numparticles > 0
            PAW = P(1:numparticles,:);
            PBK = P(numparticles+1:2*numparticles,:);
            PHARV = P(2*numparticles+1:3*numparticles,:);
            PMCD = P(3*numparticles+1:4*numparticles,:);
            PWEND = P(4*numparticles+1:5*numparticles,:);
        else
            PAW = P(1,:);
            PBK = P(2,:);
            PHARV = P(3,:);
            PMCD = P(4,:);
            PWEND = P(5,:);
        end

        % Calculate the joint probability.
        if numparticles > 0
            jointP = PAW .* PBK .* PHARV .* PMCD .* PWEND;
            pbound = 0.000000001 * ones(length(jointP),1);                                      % Establish bound.
            jointP = 1*(jointP > 0.000000001) .* jointP + 1*(jointP < 0.000000001) .* pbound;   % Use bound if underflow.

            % Average the likelihood across simulations and store into vector.
            mjointP = mean(jointP);
            lmjointP(t) = log(mjointP);

            % Resampling from Z's in time t.
            rindex = randsample(index, numparticles, true, jointP);         % Draw the index.
            lagZ = rp(rindex,:);                                            % Use the resample index to get new Z's.
        else
            lmjointP(t) = log(PAW) + log(PBK) + log(PHARV) + log(PMCD) + log(PWEND);
        end

    end

    % Store the vector of log-likelihoods over time for each market iteration.
    ll(:,market) = lmjointP;

end

% Sums the log-likelihoods across all markets and time.
sll = -sum(sum(real(ll)));

end
