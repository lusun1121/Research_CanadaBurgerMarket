% Dynamic Spillovers in the Retail Industry.                    %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.              %
% Forward simulate the structural profits for given policy.     %
% September 23, 2014.                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [APIAW, APIBK, APIHARV, APIMCD, APIWEND] = forwardsim(param, perturb, rhoX, sigmaX, data, stdX, specification, numpaths, numperiods, discount)

% Settings for simulation.
years = 36;                         % Number of years.
nummarkets = length(data)/years;    % Number of markets.
players = 5;                        % Number of players.
nexogeq = 4;                        % Number of equations in exogenous SUR process.
numcolumns = 27;                    % Number of linear coefficients to store.

% Parse first stage parameters.
[theta, mu, delta, beta, eta] = getparam1(param, specification);

% Perturb policy parameters.
theta = theta + perturb;

% Matrix to store the average of simulated profit streams for a given policy.
APIAW = zeros(nummarkets, numcolumns);
APIBK = zeros(nummarkets, numcolumns);
APIHARV = zeros(nummarkets, numcolumns);
APIMCD = zeros(nummarkets, numcolumns);
APIWEND = zeros(nummarkets, numcolumns);

% Useful vectors of ones and zeros
op = ones(numpaths, 1);
oc = ones(1, numcolumns);
zp = zeros(numpaths, 1);
zM = zeros(numpaths, nexogeq+1);

% Loop over all markets.
for market = 1:nummarkets

    % Simulation of paths for each market.
    randepsilons = normrnd(0,1,numpaths*numperiods,players);             % Draws for the epsilons in Z process.
    randu = normrnd(0,1,numpaths*numperiods,players);
    randexog = mvnrnd(zeros(nexogeq,1),sigmaX,numpaths*numperiods);       % Generates the distribution of multivariate shocks according to SUR correlation matrix.

    % Parse out the data based on market index.
    dataN = data((market-1)*years+1:market*years,:);

    % Parse out the data.
    time = dataN(:,1);           % Year.
    numAW = dataN(:,2);          % Number of A & W outlets in city.
    numBK = dataN(:,3);          % Number of Burger King outlets in city.
    numHARV = dataN(:,4);        % Number of Harvey's outlets in city.
    numMCD = dataN(:,5);         % Number of McDonald's outlets in city.
    numWEND = dataN(:,6);        % Number of Wendy's outlets in city.
    population = dataN(:,12);    % Population.
    income = dataN(:,13);        % Income.
    value = dataN(:,14);         % Property value.
    cityindex = dataN(:,20);     % Index for city.

    % Exogenous shock to city's exposure.
    greycup = dataN(:,45);              % City hosting CFL Grey Cup.
    regulation = dataN(:,46);           % Introduction of smoking regulation in city.
    minwage = dataN(:,47);              % Minimum wage.

    % Initialize X.
    D = [population(1), income(1), value(1), minwage(1)];
    Xo = kron(D,op);
    Cind = kron(cityindex(1),op);

    % Initialize storage for simulated X, Z, N, and action paths
    X = zeros(numpaths*numperiods,nexogeq);
    randparticles = zeros(numpaths*numperiods,players);
    lagnum = zeros(numpaths*numperiods,players);
    lagact = zeros(numpaths*numperiods,players);

    % Initialize storage for simulated payoff coefficients
    PIAW = zeros(numpaths,numcolumns);
    PIBK = zeros(numpaths,numcolumns);
    PIHARV = zeros(numpaths,numcolumns);
    PIMCD = zeros(numpaths,numcolumns);
    PIWEND = zeros(numpaths,numcolumns);

    % Simulate numperiods years of data
    for t = 1:numperiods

        % Parse out the draws and relevant data.
        if t == 1
            % Simulation draws
            re = randepsilons(1:numpaths,:);
            ru = randu(1:numpaths,:);

            % Inital values of lagged variables
            lagN = zeros(numpaths,players);

            % Initial values of X
            Xhatrshape = Xo;
            host = kron(greycup(1),op);
            smoke = kron(regulation(1),op);

            % Initial draws of Z
            rp = re;
            randparticles(1:numpaths,:) = rp;
        else
            % Simulation draws
            re = randepsilons((t-1)*numpaths+1:t*numpaths,:);
            rx = randexog((t-1)*numpaths+1:t*numpaths,:);
            ru = randu((t-1)*numpaths+1:t*numpaths,:);

            % Transition Z
            lagN = lagnum((t-1)*numpaths+1:t*numpaths,:);
            lagA = lagact((t-1)*numpaths+1:t*numpaths,:);
            lagZ = [rp(:,1), rp(:,2), rp(:,3), rp(:,4), rp(:,5)];
            Zhat = draw_Z(lagZ, lagA, lagN, Cind, re, mu, delta, beta, eta, specification);
            rp = [Zhat(:,1), Zhat(:,2), Zhat(:,3), Zhat(:,4), Zhat(:,5)];

            % Transition X
            Xlag = [M, zM, zM, zM; zM, M, zM, zM; zM, zM, M, zM; zM, zM, zM, M];
            rxrshape = reshape(rx,numpaths*nexogeq,1);   % Reshape the random draws
            Xhat = Xlag * rhoX + rxrshape;               % Simulate X using SUR.
            Xhatrshape = reshape(Xhat,numpaths,nexogeq); % Reshape the Xhat vector
            if t <= years
                % Stop updating GC and RI after observable period
                host = kron(greycup(t),op);                       % Perfect foresight about CFL Grey Cup.
                smoke = kron(regulation(t),op);                   % Perfect foresight about municipal smoking regulations.
            end
        end

        % Store current X and Z values
        X((t-1)*numpaths+1:t*numpaths,:) = Xhatrshape;
        randparticles((t-1)*numpaths+1:t*numpaths,:) = rp;
        M = [Xhatrshape, op];

        % Rescale variables
        pop = M(:,1) / stdX(12);
        inc = M(:,2) / stdX(13);
        val = M(:,3) / stdX(14);
        minw = M(:,4) / stdX(47);

        % Construct data matrices
        D_ALL = [ pop, inc, val, host, smoke, minw, ...
                  pop .* pop, pop .* inc, pop .* val, pop .* minw, ...
                  inc .* inc, inc .* val, inc .* minw, ...
                  val .* val, val .* minw, ...
                  minw .* minw ];
        D_interact = [ op, pop, inc, val, minw ];
        D_AW = [ D_ALL, diag(rp(:,1)) * D_interact, diag(mean(rp(:,[2,3,4,5]), 2)) * D_interact, lagN(:,1) ];
        D_BK = [ D_ALL, diag(rp(:,2)) * D_interact, diag(mean(rp(:,[1,3,4,5]), 2)) * D_interact, lagN(:,2) ];
        D_HARV = [ D_ALL, diag(rp(:,3)) * D_interact, diag(mean(rp(:,[1,2,4,5]), 2)) * D_interact, lagN(:,3) ];
        D_MCD = [ D_ALL, diag(rp(:,4)) * D_interact, diag(mean(rp(:,[1,2,3,5]), 2)) * D_interact, lagN(:,4) ];
        D_WEND = [ D_ALL, diag(rp(:,5)) * D_interact, diag(mean(rp(:,[1,2,3,4]), 2)) * D_interact, lagN(:,5) ];
        u = [ ru(:,1); ru(:,2); ru(:,3); ru(:,4); ru(:,5); ];
        D0 = zeros(size(D_AW,1), size(D_AW,2));
        D = [ D_AW D0 D0 D0 D0; D0 D_BK D0 D0 D0; D0 D0 D_HARV D0 D0; D0 D0 D0 D_MCD D0; D0 D0 D0 D0 D_WEND ];        % Compile the vectors of X and Z that affect decision.

        % Simulate actions
        [A, Xb] = forwardaction(theta,D,u);         % Simulated actions.
        A = reshape(A,numpaths,players);            % Reshape simulated action.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate the one-shot payoff.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        lagNAW = lagN(:,1);     % Parse out lagged A & W store numbers.
        lagNBK = lagN(:,2);     % Parse out lagged Burger King store numbers.
        lagNHARV = lagN(:,3);   % Parse out lagged Harvey's store numbers.
        lagNMCD = lagN(:,4);    % Parse out lagged McDonald's store numbers.
        lagNWEND = lagN(:,5);   % Parse out lagged Wendy's store numbers.

        AAW = A(:,1);           % Parse out simulated A & W action.
        ABK = A(:,2);           % Parse out simulated Burger King action.
        AHARV = A(:,3);         % Parse out simulated Harvey's action.
        AMCD = A(:,4);          % Parse out simulated McDonald's action.
        AWEND = A(:,5);         % Parse out simulated Wendy's action.

        NAW = lagNAW + AAW;       % Parse out simulated A & W store number.
        NBK = lagNBK + ABK;       % Parse out simulated Burger King store number.
        NHARV = lagNHARV + AHARV; % Parse out simulated Harvey's store number.
        NMCD = lagNMCD + AMCD;    % Parse out simulated McDonald's store number.
        NWEND = lagNWEND + AWEND; % Parse out simulated Wendy's store number.

        % Non-negativity constraint for store numbers.
        NAW = 1*(NAW > 0).*NAW;
        NBK = 1*(NBK > 0).*NBK;
        NHARV = 1*(NHARV > 0).*NHARV;
        NMCD = 1*(NMCD > 0).*NMCD;
        NWEND = 1*(NWEND > 0).*NWEND;

        NAWSQ = (NAW ./ 100).^2;
        NBKSQ = (NBK ./ 100).^2;
        NHARVSQ = (NHARV ./ 100).^2;
        NMCDSQ = (NMCD ./ 100).^2;
        NWENDSQ = (NWEND ./ 100).^2;

        activeAW = 1*(AAW > 0 | lagNAW > 0);
        activeBK = 1*(ABK > 0 | lagNBK > 0);
        activeHARV = 1*(AHARV > 0 | lagNHARV > 0);
        activeMCD = 1*(AMCD > 0 | lagNMCD > 0);
        activeWEND = 1*(AWEND > 0 | lagNWEND > 0);

        logNAW = log(NAW+1); logNAW(~activeAW) = 0;
        logNBK = log(NBK+1); logNBK(~activeBK) = 0;
        logNHARV = log(NHARV+1); logNHARV(~activeHARV) = 0;
        logNMCD = log(NMCD+1); logNMCD(~activeMCD) = 0;
        logNWEND = log(NWEND+1); logNWEND(~activeWEND) = 0;

        NcompAW = NBK + NHARV + NMCD + NWEND;   % Competitors of A & W.
        NcompBK = NAW + NHARV + NMCD + NWEND;   % Competitors of Burger King.
        NcompHARV = NAW + NBK + NMCD + NWEND;   % Competitors of Harvey's.
        NcompMCD = NAW + NBK + NHARV + NWEND;   % Competitors of McDonald's.
        NcompWEND = NAW + NBK + NHARV + NMCD;   % Competitors of Wendy's.

        activecompAW = 1*(NcompAW > 0);
        activecompBK = 1*(NcompBK > 0);
        activecompHARV = 1*(NcompHARV > 0);
        activecompMCD = 1*(NcompMCD > 0);
        activecompWEND = 1*(NcompWEND > 0);

        logNcompAW = log(NcompAW+1); logNcompAW(~activecompAW) = 0;
        logNcompBK = log(NcompBK+1); logNcompBK(~activecompBK) = 0;
        logNcompHARV = log(NcompHARV+1);logNcompHARV(~activecompHARV) = 0;
        logNcompMCD = log(NcompMCD+1); logNcompMCD(~activecompMCD) = 0;
        logNcompWEND = log(NcompWEND+1); logNcompWEND(~activecompWEND) = 0;

        enterAW = 1*(lagNAW == 0 & AAW > 0);         % Entry term for A & W.
        enterBK = 1*(lagNBK == 0 & ABK > 0);         % Entry term for Burger King.
        enterHARV = 1*(lagNHARV == 0 & AHARV > 0);   % Entry term for Harvey's.
        enterMCD = 1*(lagNMCD == 0 & AMCD > 0);      % Entry term for McDonald's.
        enterWEND = 1*(lagNWEND == 0 & AWEND > 0);   % Entry term for Wendy's.

        expandAW = 1*(AAW > 0) .* AAW;         % Expansion term for A & W.
        expandBK = 1*(ABK > 0) .* ABK;         % Expansion term for Burger King.
        expandHARV = 1*(AHARV > 0) .* AHARV;   % Expansion term for Harvey's.
        expandMCD = 1*(AMCD > 0) .* AMCD;      % Expansion term for McDonald's.
        expandWEND = 1*(AWEND > 0) .* AWEND;   % Expansion term for Wendy's.

        contractAW = 1*(lagNAW > 0 & AAW < 0) .* -AAW;         % Contraction term for A & W.
        contractBK = 1*(lagNAW > 0 & ABK < 0) .* -ABK;         % Contraction term for Burger King.
        contractHARV = 1*(lagNAW > 0 & AHARV < 0) .* -AHARV;   % Contraction term for Harvey's.
        contractMCD = 1*(lagNAW > 0 & AMCD < 0) .* -AMCD;      % Contraction term for McDonald's.
        contractWEND = 1*(lagNAW > 0 & AWEND < 0) .* -AWEND;   % Contraction term for Wendy's.

        % Calculate the variable profits for each chain.
        simVAW = [ pop, inc, val, host, smoke, minw, ...
                   pop .* NAW, inc .* NAW, val .* NAW, host .* NAW, smoke .* NAW, minw .* NAW, ...
                   rp(:,1), NcompAW, logNcompAW, zp, NBK, NHARV, NMCD, NWEND, ...
                   NAW, NAWSQ, logNAW, ...
                   enterAW, expandAW, contractAW, ru(:,1) ];
        simVBK = [ pop, inc, val, host, smoke, minw, ...
                   pop .* NBK, inc .* NBK, val .* NBK, host .* NBK, smoke .* NBK, minw .* NBK, ...
                   rp(:,2), NcompBK, logNcompBK, NAW, zp, NHARV, NMCD, NWEND, ...
                   NBK, NBKSQ, logNBK, ...
                   enterBK, expandBK, contractBK, ru(:,2) ];
        simVHARV = [ pop, inc, val, host, smoke, minw, ...
                     pop .* NHARV, inc .* NHARV, val .* NHARV, host .* NHARV, smoke .* NHARV, minw .* NHARV, ...
                     rp(:,3), NcompHARV, logNcompHARV, NAW, NBK, zp, NMCD, NWEND, ...
                     NHARV, NHARVSQ, logNHARV, ...
                     enterHARV, expandHARV, contractHARV, ru(:,3) ];
        simVMCD = [ pop, inc, val, host, smoke, minw, ...
                    pop .* NMCD, inc .* NMCD, val .* NMCD, host .* NMCD, smoke .* NMCD, minw .* NMCD, ...
                    rp(:,4), NcompMCD, logNcompMCD, NAW, NBK, NHARV, zp, NWEND, ...
                    NMCD, NMCDSQ, logNMCD, ...
                    enterMCD, expandMCD, contractMCD, ru(:,4) ];
        simVWEND = [ pop, inc, val, host, smoke, minw, ...
                     pop .* NWEND, inc .* NWEND, val .* NWEND, host .* NWEND, smoke .* NWEND, minw .* NWEND, ...
                     rp(:,5), NcompWEND, logNcompWEND, NAW, NBK, NHARV, NMCD, zp, ...
                     NWEND, NWENDSQ, logNWEND, ...
                     enterWEND, expandWEND, contractWEND, ru(:,5) ];

        actAW = kron(activeAW, oc);
        actBK = kron(activeBK, oc);
        actHARV = kron(activeHARV, oc);
        actMCD = kron(activeMCD, oc);
        actWEND = kron(activeWEND, oc);

        simVAW = simVAW .* actAW;
        simVBK = simVBK .* actBK;
        simVHARV = simVHARV .* actHARV;
        simVMCD = simVMCD .* actMCD;
        simVWEND = simVWEND .* actWEND;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Update the discounted sum of profits with the current net profits.
        disc = discount^(t-1);
        PIAW = PIAW + disc*simVAW;
        PIBK = PIBK + disc*simVBK;
        PIHARV = PIHARV + disc*simVHARV;
        PIMCD = PIMCD + disc*simVMCD;
        PIWEND = PIWEND + disc*simVWEND;

        if t < numperiods
            % Store total number of stores for next period.
            lagnum(t*numpaths+1:(t+1)*numpaths,:) = [ NAW, NBK, NHARV, NMCD, NWEND ];

            % Store net number of stores added next period.
            lagact(t*numpaths+1:(t+1)*numpaths,:) = [ AAW, ABK, AHARV, AMCD, AWEND ];
        end
    end

    % Average the discounted profit streams across all simulations.
    APIAW(market,:) = mean(PIAW);
    APIBK(market,:) = mean(PIBK);
    APIHARV(market,:) = mean(PIHARV);
    APIMCD(market,:) = mean(PIMCD);
    APIWEND(market,:) = mean(PIWEND);

end

end
