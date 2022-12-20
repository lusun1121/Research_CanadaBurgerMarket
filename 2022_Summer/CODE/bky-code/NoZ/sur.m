% Dynamic Spillovers in the Retail Industry.        %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.  %
% SUR model estimation.                             %
% November 25, 2012.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rhoX, sigmaX, se] = sur(data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data set-up step.                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse out the data.
time = data(:,1);           % Year.
numAW = data(:,2);          % Number of A & W outlets in city.
numBK = data(:,3);          % Number of Burger King outlets in city.
numHARV = data(:,4);        % Number of Harvey's outlets in city.
numMCD = data(:,5);         % Number of McDonald's outlets in city.
numWEND = data(:,6);        % Number of Wendy's outlets in city.
lagnumAW = data(:,7);       % Lagged number of A & W outlets in city.
lagnumBK = data(:,8);       % Lagged number of Burger King outlets in city.
lagnumHARV = data(:,9);     % Lagged number of Harvey's outlets in city.
lagnumMCD = data(:,10);     % Lagged number of McDonald's outlets in city.
lagnumWEND = data(:,11);    % Lagged number of Wendy's outlets in city.
population = data(:,12);    % Population.
income = data(:,13);        % Income.
value = data(:,14);         % Property value.
aAW = data(:,15);           % Change in A & W outlets.
aBK = data(:,16);           % Change in Burger King outlets.
aHARV = data(:,17);         % Change in Harvey's outlets.
aMCD = data(:,18);          % Change in McDonald's outlets.
aWEND = data(:,19);         % Change in Wendy's outlets.
cityindex = data(:,20);     % Index for city.
timeindex = data(:,21);     % Index for time.
lagnum2AW = data(:,22);     % Squared lagged number of A & W outlets in city.
lagnum2BK = data(:,23);     % Squared lagged number of Burger King outlets in city.
lagnum2HARV = data(:,24);   % Squared lagged number of Harvey's outlets in city.
lagnum2MCD = data(:,25);    % Squared lagged number of McDonald's outlets in city.
lagnum2WEND = data(:,26);   % Squared lagged number of Wendy's outlets in city.
lagtotalAW = data(:,27);    % National level of A & W stores.
lagtotalBK = data(:,28);    % National level of Burger King stores.
lagtotalHARV = data(:,29);  % National level of Harvey's stores.
lagtotalMCD = data(:,30);   % National level of McDonald's stores.
lagtotalWEND = data(:,31);  % National level of Wendy's stores.
lagtotal2AW = data(:,32);   % Squared national level of A & W stores.
lagtotal2BK = data(:,33);   % Squared national level of Burger King stores.
lagtotal2HARV = data(:,34); % Squared national level of Harvey's stores.
lagtotal2MCD = data(:,35);  % Squared national level of McDonald's stores.
lagtotal2WEND = data(:,36); % Squared national level of Wendy's stores.

% Variables used for VAR regression.
lagpopulation = data(:,37);
lagincome = data(:,38);
lagvalue = data(:,39);
totalAW = data(:,40);
totalBK = data(:,41);
totalHARV = data(:,42);
totalMCD = data(:,43);
totalWEND = data(:,44);
minwage = data(:,47);
lagminwage = data(:,48);

% Create new state variables.
lagcityotherAW = lagnumBK + lagnumHARV + lagnumMCD + lagnumWEND;    % Lagged number of A & W's competitors in same city.
lagcityotherBK = lagnumAW + lagnumHARV + lagnumMCD + lagnumWEND;    % Lagged number of Burger King's competitors in same city.
lagcityotherHARV = lagnumAW + lagnumBK + lagnumMCD + lagnumWEND;    % Lagged number of Harvey's competitors in same city.
lagcityotherMCD = lagnumAW + lagnumBK + lagnumHARV + lagnumWEND;    % Lagged number of McDonald's competitors in same city.
lagcityotherWEND = lagnumAW + lagnumBK + lagnumHARV + lagnumMCD;    % Lagged number of Wendy's competitors in same city.
lagcityother2AW = lagcityotherAW.^2;                                % Squared number of A & W's competitors in same city.
lagcityother2BK = lagcityotherBK.^2;                                % Squared number of Burger King's competitors in same city.
lagcityother2HARV = lagcityotherHARV.^2;                            % Squared number of Harvey's competitors in same city.
lagcityother2MCD = lagcityotherMCD.^2;                              % Squared number of McDonald's competitors in same city.
lagcityother2WEND = lagcityotherWEND.^2;                            % Squared number of Wendy's competitors in same city.
lagtotalotherAW = lagtotalBK + lagtotalHARV + lagtotalMCD + lagtotalWEND;   % Lagged number of A & W's competitors in country.
lagtotalotherBK = lagtotalAW + lagtotalHARV + lagtotalMCD + lagtotalWEND;   % Lagged number of Burger King's competitors in country.
lagtotalotherHARV = lagtotalAW + lagtotalBK + lagtotalMCD + lagtotalWEND;   % Lagged number of Harvey's competitors in country.
lagtotalotherMCD = lagtotalAW + lagtotalBK + lagtotalHARV + lagtotalWEND;   % Lagged number of McDonald's competitors in country.
lagtotalotherWEND = lagtotalAW + lagtotalBK + lagtotalHARV + lagtotalMCD;   % Lagged number of Wendy's competitors in country.
lagtotalother2AW = lagtotalotherAW.^2;                                      % Squared number of A & W's competitors in country.
lagtotalother2BK = lagtotalotherBK.^2;                                      % Squared number of Burger King's competitors in country.
lagtotalother2HARV = lagtotalotherHARV.^2;                                  % Squared number of Harvey's competitors in country.
lagtotalother2MCD = lagtotalotherMCD.^2;                                    % Squared number of McDonald's competitors in country.
lagtotalother2WEND = lagtotalotherWEND.^2;                                  % Squared number of Wendy's competitors in country.
mt = size(data,1);                                                          % Total number of observations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reshape variables to remove first period for each market
years = 36;                         % Number of years.
markets = length(data) / years;     % Number of markets.

popr = reshape(population, years, markets);
incr = reshape(income, years, markets);
valr = reshape(value, years, markets);
minwr = reshape(minwage, years, markets);
lagpopr = reshape(lagpopulation, years, markets);
lagincr = reshape(lagincome, years, markets);
lagvalr = reshape(lagvalue, years, markets);
lagminwr = reshape(lagminwage, years, markets);

popr = popr(2:end,:);
incr = incr(2:end,:);
valr = valr(2:end,:);
minwr = minwr(2:end,:);
lagpopr = lagpopr(2:end,:);
lagincr = lagincr(2:end,:);
lagvalr = lagvalr(2:end,:);
lagminwr = lagminwr(2:end,:);

pop = reshape(popr, markets*(years-1), 1);
inc = reshape(incr, markets*(years-1), 1);
val = reshape(valr, markets*(years-1), 1);
minw = reshape(minwr, markets*(years-1), 1);
lagpop = reshape(lagpopr, markets*(years-1), 1);
laginc = reshape(lagincr, markets*(years-1), 1);
lagval = reshape(lagvalr, markets*(years-1), 1);
lagminw = reshape(lagminwr, markets*(years-1), 1);

% Specify the exogenous variables for SUR estimation.
Y = [pop, inc, val, minw];
lagY = [lagpop, laginc, lagval, lagminw];

[results, sigmaX] = sureg(Y,lagY);

rhoX = results(:,1);
se = results(:,2);

end
