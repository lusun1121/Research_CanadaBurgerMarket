% Dynamic Spillovers in the Retail Industry.        %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.  %
% Draw from the distribution of Z.                  %
% March 7, 2013.                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Zhat = draw_Z(lagZ, lagA, lagN, cityindex, epsilon, mu, delta, beta, eta, specification)

% Dimensions
numpaths = size(lagZ, 1);
o = ones(numpaths, 1);

% Parse out the lagged variables.
lagZAW = lagZ(:,1);
lagZBK = lagZ(:,2);
lagZHARV = lagZ(:,3);
lagZMCD = lagZ(:,4);
lagZWEND = lagZ(:,5);

lagAAW = lagA(:,1);
lagABK = lagA(:,2);
lagAHARV = lagA(:,3);
lagAMCD = lagA(:,4);
lagAWEND = lagA(:,5);

lagNAW = lagN(:,1);
lagNBK = lagN(:,2);
lagNHARV = lagN(:,3);
lagNMCD = lagN(:,4);
lagNWEND = lagN(:,5);

% Parse out the iid draws.
epsilonAW = epsilon(:,1);
epsilonBK = epsilon(:,2);
epsilonHARV = epsilon(:,3);
epsilonMCD = epsilon(:,4);
epsilonWEND = epsilon(:,5);

% City fixed effects.
ce = eta(cityindex);

% Parse out firm-specific intercepts.
muAW = mu(1) * o;
muBK = mu(2) * o;
muHARV = mu(3) * o;
muMCD = mu(4) * o;
muWEND = mu(5) * o;

% Parse out firm-specific deltas and drifts.
deltaAW = delta(1);
deltaBK = delta(2);
deltaHARV = delta(3);
deltaMCD = delta(4);
deltaWEND = delta(5);

% Parse out the firm-specific betas.
beta1AW = beta(1);
beta2AW = beta(2);
beta3AW = beta(3);
beta4AW = beta(4);
beta5AW = beta(5);
sigmaAW = beta(6);

beta1BK = beta(7);
beta2BK = beta(8);
beta3BK = beta(9);
beta4BK = beta(10);
beta5BK = beta(11);
sigmaBK = beta(12);

beta1HARV = beta(13);
beta2HARV = beta(14);
beta3HARV = beta(15);
beta4HARV = beta(16);
beta5HARV = beta(17);
sigmaHARV = beta(18);

beta1MCD = beta(19);
beta2MCD = beta(20);
beta3MCD = beta(21);
beta4MCD = beta(22);
beta5MCD = beta(23);
sigmaMCD = beta(24);

beta1WEND = beta(25);
beta2WEND = beta(26);
beta3WEND = beta(27);
beta4WEND = beta(28);
beta5WEND = beta(29);
sigmaWEND = beta(30);

% Generate the new Z's.
spilloverAW = (1*(lagAAW > 0)).*(beta1AW * lagAAW) + (1*(lagAAW < 0)).*(beta2AW * lagAAW) + beta3AW * lagNAW + beta4AW * lagNAW.^2 + beta5AW * lagNAW.^3;
spilloverBK = (1*(lagABK > 0)).*(beta1BK * lagABK) + (1*(lagABK < 0)).*(beta2BK * lagABK) + beta3BK * lagNBK + beta4BK * lagNBK.^2 + beta5BK * lagNBK.^3;
spilloverHARV = (1*(lagAHARV > 0)).*(beta1HARV * lagAHARV) + (1*(lagAHARV < 0)).*(beta2HARV * lagAHARV) + beta3HARV * lagNHARV + beta4HARV * lagNHARV.^2 + beta5HARV * lagNHARV.^3;
spilloverMCD = (1*(lagAMCD > 0)).*(beta1MCD * lagAMCD) + (1*(lagAMCD < 0)).*(beta2MCD * lagAMCD) + beta3MCD * lagNMCD + beta4MCD * lagNMCD.^2 + beta5MCD * lagNMCD.^3;
spilloverWEND = (1*(lagAWEND > 0)).*(beta1WEND * lagAWEND) + (1*(lagAWEND < 0)).*(beta2WEND * lagAWEND) + beta3WEND * lagNWEND + beta4WEND * lagNWEND.^2 + beta5WEND * lagNWEND.^3;

ZAW = muAW + deltaAW * lagZAW + spilloverAW + ce + sigmaAW * epsilonAW;
ZBK = muBK + deltaBK * lagZBK + spilloverBK +ce + sigmaBK * epsilonBK;
ZHARV = muHARV + deltaHARV * lagZHARV + spilloverHARV + ce + sigmaHARV * epsilonHARV;
ZMCD = muMCD + deltaMCD * lagZMCD + spilloverMCD + ce + sigmaMCD * epsilonMCD;
ZWEND = muWEND + deltaWEND * lagZWEND + spilloverWEND + ce + sigmaWEND * epsilonWEND;

Zhat = [ZAW, ZBK, ZHARV, ZMCD, ZWEND];

end
