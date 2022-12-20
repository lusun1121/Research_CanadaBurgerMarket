% Dynamic Spillovers in the Retail Industry.                    %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.              %
% Forward simulate the structural profits for given policy.     %
% January 6, 2015                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N_mean, N_bias, N_mse, Zstore, PIstore] = forwardsiminfo(param, rhoX, sigmaX, data, stdX, specification, numpaths, numperiods, discount, paramS, secondspec, counterfactual)

% Settings for simulation.
years = 36;                         % Number of years.
nummarkets = length(data)/years;    % Number of markets.
players = 5;                        % Number of players.
nexogeq = 4;                        % Number of equations in exogenous SUR process.
numcolumns = 27;                    % Number of linear coefficients to store.

% Parse first stage parameters.
[theta, mu, delta, beta, eta] = getparam1(param, specification);

% Parse second stage parameters.
thetaS = getparam2(paramS, secondspec);

% Matrix to store the simulated N's and PI's.
N_mean = zeros(nummarkets*years,players);
N_var = zeros(nummarkets*years,players);
N_bias = zeros(nummarkets*years,players);
N_mse = zeros(nummarkets*years,players);
PIstore = zeros(nummarkets*years,players);
Zstore = zeros(nummarkets*years,players);

% Useful vectors of ones and zeros
op = ones(numpaths, 1);
oc = ones(1, numcolumns);
zp = zeros(numpaths, 1);
zM = zeros(numpaths, nexogeq+1);

% Loop over all markets.
for market = 1:nummarkets

    % Simulation of paths for each market.
    stationarymeanAW = (mu(1)+eta(market))/(1-delta(1));
    stationarymeanBK = (mu(2)+eta(market))/(1-delta(2));
    stationarymeanHARV = (mu(3)+eta(market))/(1-delta(3));
    stationarymeanMCD = (mu(4)+eta(market))/(1-delta(4));
    stationarymeanWEND = (mu(5)+eta(market))/(1-delta(5));
    stationaryvarAW = (beta(6)^2)/(1-delta(1)^2);
    stationaryvarBK = (beta(12)^2)/(1-delta(2)^2);
    stationaryvarHARV = (beta(18)^2)/(1-delta(3)^2);
    stationaryvarMCD = (beta(24)^2)/(1-delta(4)^2);
    stationaryvarWEND = (beta(30)^2)/(1-delta(5)^2);

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
    %lagnumAW = dataN(:,7);       % Lagged number of A & W outlets in city.
    %lagnumBK = dataN(:,8);       % Lagged number of Burger King outlets in city.
    %lagnumHARV = dataN(:,9);     % Lagged number of Harvey's outlets in city.
    %lagnumMCD = dataN(:,10);     % Lagged number of McDonald's outlets in city.
    %lagnumWEND = dataN(:,11);    % Lagged number of Wendy's outlets in city.
    population = dataN(:,12);    % Population.
    income = dataN(:,13);        % Income.
    value = dataN(:,14);         % Property value.
    %aAW = dataN(:,15);           % Change in A & W outlets.
    %aBK = dataN(:,16);           % Change in Burger King outlets.
    %aHARV = dataN(:,17);         % Change in Harvey's outlets.
    %aMCD = dataN(:,18);          % Change in McDonald's outlets.
    %aWEND = dataN(:,19);         % Change in Wendy's outlets.
    cityindex = dataN(:,20);     % Index for city.
    %timeindex = dataN(:,21);     % Index for time.
    %lagnum2AW = dataN(:,22);     % Squared lagged number of A & W outlets in city.
    %lagnum2BK = dataN(:,23);     % Squared lagged number of Burger King outlets in city.
    %lagnum2HARV = dataN(:,24);   % Squared lagged number of Harvey's outlets in city.
    %lagnum2MCD = dataN(:,25);    % Squared lagged number of McDonald's outlets in city.
    %lagnum2WEND = dataN(:,26);   % Squared lagged number of Wendy's outlets in city.
    %lagtotalAW = dataN(:,27);    % National level of A & W stores.
    %lagtotalBK = dataN(:,28);    % National level of Burger King stores.
    %lagtotalHARV = dataN(:,29);  % National level of Harvey's stores.
    %lagtotalMCD = dataN(:,30);   % National level of McDonald's stores.
    %lagtotalWEND = dataN(:,31);  % National level of Wendy's stores.
    %lagtotal2AW = dataN(:,32);   % Squared national level of A & W stores.
    %lagtotal2BK = dataN(:,33);   % Squared national level of Burger King stores.
    %lagtotal2HARV = dataN(:,34); % Squared national level of Harvey's stores.
    %lagtotal2MCD = dataN(:,35);  % Squared national level of McDonald's stores.
    %lagtotal2WEND = dataN(:,36); % Squared national level of Wendy's stores.

    % Variables used for SUR regression.
    %lagpopulation = dataN(:,37);
    %lagincome = dataN(:,38);
    %lagvalue = dataN(:,39);
    %totalAW = dataN(:,40);
    %totalBK = dataN(:,41);
    %totalHARV = dataN(:,42);
    %totalMCD = dataN(:,43);
    %totalWEND = dataN(:,44);

    % Exogenous shock to city's exposure.
    greycup = dataN(:,45);              % City hosting CFL Grey Cup.
    regulation = dataN(:,46);           % Introduction of smoking regulation in city.
    minwage = dataN(:,47);              % Minimum wage.
    %lagminwage = dataN(:,48);           % Lagged minimum wage.

    % Lagged actions.
    %lagaAW = dataN(:,49);
    %lagaBK = dataN(:,50);
    %lagaHARV = dataN(:,51);
    %lagaMCD = dataN(:,52);
    %lagaWEND = dataN(:,53);

    % Create new state variables.
    %lagcityotherAW = lagnumBK + lagnumHARV + lagnumMCD + lagnumWEND;    % Lagged number of A & W's competitors in same city.
    %lagcityotherBK = lagnumAW + lagnumHARV + lagnumMCD + lagnumWEND;    % Lagged number of Burger King's competitors in same city.
    %lagcityotherHARV = lagnumAW + lagnumBK + lagnumMCD + lagnumWEND;    % Lagged number of Harvey's competitors in same city.
    %lagcityotherMCD = lagnumAW + lagnumBK + lagnumHARV + lagnumWEND;    % Lagged number of McDonald's competitors in same city.
    %lagcityotherWEND = lagnumAW + lagnumBK + lagnumHARV + lagnumMCD;    % Lagged number of Wendy's competitors in same city.
    %lagcityother2AW = lagcityotherAW.^2;                                % Squared number of A & W's competitors in same city.
    %lagcityother2BK = lagcityotherBK.^2;                                % Squared number of Burger King's competitors in same city.
    %lagcityother2HARV = lagcityotherHARV.^2;                            % Squared number of Harvey's competitors in same city.
    %lagcityother2MCD = lagcityotherMCD.^2;                              % Squared number of McDonald's competitors in same city.
    %lagcityother2WEND = lagcityotherWEND.^2;                            % Squared number of Wendy's competitors in same city.
    %lagtotalotherAW = lagtotalBK + lagtotalHARV + lagtotalMCD + lagtotalWEND;   % Lagged number of A & W's competitors in country.
    %lagtotalotherBK = lagtotalAW + lagtotalHARV + lagtotalMCD + lagtotalWEND;   % Lagged number of Burger King's competitors in country.
    %lagtotalotherHARV = lagtotalAW + lagtotalBK + lagtotalMCD + lagtotalWEND;   % Lagged number of Harvey's competitors in country.
    %lagtotalotherMCD = lagtotalAW + lagtotalBK + lagtotalHARV + lagtotalWEND;   % Lagged number of McDonald's competitors in country.
    %lagtotalotherWEND = lagtotalAW + lagtotalBK + lagtotalHARV + lagtotalMCD;   % Lagged number of Wendy's competitors in country.
    %lagtotalother2AW = lagtotalotherAW.^2;                                      % Squared number of A & W's competitors in country.
    %lagtotalother2BK = lagtotalotherBK.^2;                                      % Squared number of Burger King's competitors in country.
    %lagtotalother2HARV = lagtotalotherHARV.^2;                                  % Squared number of Harvey's competitors in country.
    %lagtotalother2MCD = lagtotalotherMCD.^2;                                    % Squared number of McDonald's competitors in country.
    %lagtotalother2WEND = lagtotalotherWEND.^2;                                  % Squared number of Wendy's competitors in country.
    %mt = size(dataN,1);                                                         % Total number of observations.

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
            NinitAW = kron(numAW(1),op);
            NinitBK = kron(numBK(1),op);
            NinitHARV = kron(numHARV(1),op);
            NinitMCD = kron(numMCD(1),op);
            NinitWEND = kron(numWEND(1),op);
            lagN = [NinitAW, NinitBK, NinitHARV, NinitMCD, NinitWEND];
            %lagN = zeros(numpaths,players);

            % Initial values of X
            Xhatrshape = Xo;
            host = kron(greycup(1),op);
            smoke = kron(regulation(1),op);
            pop = kron(population(1),op) / stdX(12);
            inc = kron(income(1),op) / stdX(13);
            val = kron(value(1),op) / stdX(14);
            minw = kron(minwage(1),op) / stdX(47);
            
            % Initial draws of Z
            rpAW = stationarymeanAW + sqrt(stationaryvarAW) * re(:,1);
            rpBK = stationarymeanBK + sqrt(stationaryvarBK) * re(:,2);
            rpHARV = stationarymeanHARV + sqrt(stationaryvarHARV) * re(:,3);
            rpMCD = stationarymeanMCD + sqrt(stationaryvarMCD) * re(:,4);
            rpWEND = stationarymeanWEND + sqrt(stationaryvarWEND) * re(:,5);
            rp = [rpAW, rpBK, rpHARV, rpMCD, rpWEND];
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

            % Counterfactuals
            if t == 31

                % Counterfactuals 1-5: a positive, one-unit shock to Z in 2000
                if counterfactual > 0 && counterfactual < 6
                    rp(:,counterfactual) = rp(:,counterfactual) + ones(numpaths,1);
                end

                % Counterfactuals -5 to -1: a negative, one-unit shock to Z in 2000
                if counterfactual > -6 && counterfactual < 0
                    rp(:,-counterfactual) = rp(:,-counterfactual) - ones(numpaths,1);
                end

                %  Counterfactuals 1-5: add a one standard deviation shock in 2000
                % if counterfactual == 1
                %     rp(:,1) = rp(:,1) + sqrt(stationaryvarAW);
                % elseif counterfactual == 2
                %     rp(:,2) = rp(:,2) + sqrt(stationaryvarBK);
                % elseif counterfactual == 3
                %     rp(:,3) = rp(:,3) + sqrt(stationaryvarHARV);
                % elseif counterfactual == 4
                %     rp(:,4) = rp(:,4) + sqrt(stationaryvarMCD);
                % elseif counterfactual == 5
                %     rp(:,5) = rp(:,5) + sqrt(stationaryvarWEND);
                % end

                % Counterfactuals 11-15: remove 20% of stores in 2000
                if counterfactual == 11
                    lagN(:,1) = 0.8 * lagN(:,1);
                elseif counterfactual == 12
                    lagN(:,2) = 0.8 * lagN(:,2);
                elseif counterfactual == 13
                    lagN(:,3) = 0.8 * lagN(:,3);
                elseif counterfactual == 14
                    lagN(:,4) = 0.8 * lagN(:,4);
                elseif counterfactual == 15
                    lagN(:,5) = 0.8 * lagN(:,5);
                end

                % Counterfactuals 21-25: increase lagged stores by 20% in 2000
                if counterfactual == 21
                    lagN(:,1) = 1.2 * lagN(:,1);
                elseif counterfactual == 22
                    lagN(:,2) = 1.2 * lagN(:,2);
                elseif counterfactual == 23
                    lagN(:,3) = 1.2 * lagN(:,3);
                elseif counterfactual == 24
                    lagN(:,4) = 1.2 * lagN(:,4);
                elseif counterfactual == 25
                    lagN(:,5) = 1.2 * lagN(:,5);
                end
            end

            % Transition X
            Xlag = [M, zM, zM, zM; zM, M, zM, zM; zM, zM, M, zM; zM, zM, zM, M];
            rxrshape = reshape(rx,numpaths*nexogeq,1);   % Reshape the random draws
            Xhat = Xlag * rhoX + rxrshape;               % Simulate X using SUR.
            Xhatrshape = reshape(Xhat,numpaths,nexogeq); % Reshape the Xhat vector
            if t <= years
                % Stop updating GC and RI after observable period
                host = kron(greycup(t),op);                       % Perfect foresight about CFL Grey Cup.
                smoke = kron(regulation(t),op);                   % Perfect foresight about municipal smoking regulations.
                pop = kron(population(t),op) / stdX(12);
                inc = kron(income(t),op) / stdX(13);
                val = kron(value(t),op) / stdX(14);
                minw = kron(minwage(t),op) / stdX(47);
            end
        end

        % Store current X and Z values
        X((t-1)*numpaths+1:t*numpaths,:) = Xhatrshape;
        randparticles((t-1)*numpaths+1:t*numpaths,:) = rp;
        M = [Xhatrshape, op];

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
        [A, Xb] = forwardactioncf(theta,D,u);       % Simulated actions.
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

        thetaSr = reshape(thetaS,numcolumns,5);
        PIAWtheta = simVAW * thetaSr(:,1);
        PIBKtheta = simVBK * thetaSr(:,2);
        PIHARVtheta = simVHARV * thetaSr(:,3);
        PIMCDtheta = simVMCD * thetaSr(:,4);
        PIWENDtheta = simVWEND * thetaSr(:,5);
        
        if t < numperiods
            % Store total number of stores for next period.
            lagnum(t*numpaths+1:(t+1)*numpaths,:) = [ NAW, NBK, NHARV, NMCD, NWEND ];

            % Store net number of stores added next period.
            lagact(t*numpaths+1:(t+1)*numpaths,:) = [ AAW, ABK, AHARV, AMCD, AWEND ];
        end
    
        % Store the simulated store numbers.    
        N_mean((market-1)*years + t,:) = [mean(NAW), mean(NBK), mean(NHARV), mean(NMCD), mean(NWEND)];
        N_var((market-1)*years + t,:) = [var(NAW), var(NBK), var(NHARV), var(NMCD), var(NWEND)];
        N_bias((market-1)*years + t,:) = N_mean((market-1)*years + t,:) - [ numAW(t), numBK(t), numHARV(t), numMCD(t), numWEND(t) ];
        N_mse((market-1)*years + t,:) = N_bias((market-1)*years + t,:).^2 + N_var((market-1)*years + t,:);
        PIstore((market-1)*years + t,:) = [mean(PIAWtheta), mean(PIBKtheta), mean(PIHARVtheta), mean(PIMCDtheta), mean(PIWENDtheta)];
        Zstore((market-1)*years + t,:) = [mean(rp(:,1)), mean(rp(:,2)), mean(rp(:,3)), mean(rp(:,4)), mean(rp(:,5))];
    end

end

end
