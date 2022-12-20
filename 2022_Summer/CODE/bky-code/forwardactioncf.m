% Dynamic Spillovers in the Retail Industry.        %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.  %
% Simulates the choices for each simulation.        %
% November 14, 2012.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, Xb] = forwardactioncf(theta,X,u)

grid = [-5 -2 -1 0 1 2 5];
ng = length(grid);

% Chain-specific means
outcomes = [  -7.00000,  -7.00000, -14.00000,  -7.00000,  -7.00000;
              -2.28571,  -2.00000,  -2.37500,  -2.50000,  -2.66667;
              -1.00000,  -1.00000,  -1.00000,  -1.00000,  -1.00000;
               0.00000,   0.00000,   0.00000,   0.00000,   0.00000;
               1.00000,   1.00000,   1.00000,   1.00000,   1.00000;
               2.00000,   2.00000,   2.00000,   2.00000,   2.00000;
               3.48000,   3.58333,   3.40000,   3.64615,   3.25000;
               9.00000,   7.33333,   9.00000,  10.35714,   6.00000; ];

% Dimensions
numplayers = 5;                         % Number of firms.
nx = size(X,2);                         % Number of state variables.
numpar = length(theta);                 % Number of parameters.
nbeta = numpar / numplayers - ng;       % Number of beta parameters per firm.

% Parsing out the parameter vector.
betaAW = theta(1:nbeta);
cutoffAW = theta(nbeta+1:nbeta+ng);
betaBK = theta(nbeta+ng+1:2*nbeta+ng);
cutoffBK = theta(2*nbeta+ng+1:2*nbeta+2*ng);
betaHARV = theta(2*nbeta+2*ng+1:3*nbeta+2*ng);
cutoffHARV = theta(3*nbeta+2*ng+1:3*nbeta+3*ng);
betaMCD = theta(3*nbeta+3*ng+1:4*nbeta+3*ng);
cutoffMCD = theta(4*nbeta+3*ng+1:4*nbeta+4*ng);
betaWEND = theta(4*nbeta+4*ng+1:5*nbeta+4*ng);
cutoffWEND = theta(5*nbeta+4*ng+1:5*nbeta+5*ng);
beta = [ betaAW; betaBK; betaHARV; betaMCD; betaWEND; ];

% Calculating contribution of state variables to the latent payoff.
Xb = X*beta + u;
N = length(Xb)/5;
XbAW = Xb(1:N,:);
XbBK = Xb(N+1:2*N,:);
XbHARV = Xb(2*N+1:3*N,:);
XbMCD = Xb(3*N+1:4*N,:);
XbWEND = Xb(4*N+1:5*N,:);

% Initialize matrix of actions.
AAW = zeros(N,ng);
ABK = zeros(N,ng);
AHARV = zeros(N,ng);
AMCD = zeros(N,ng);
AWEND = zeros(N,ng);

% Generate the simulated expansion/contraction actions.
for i = 1:ng+1
    if i == 1
        AAW(:,i) = 1*(XbAW <= cutoffAW(1)) * outcomes(1,1);
        ABK(:,i) = 1*(XbBK <= cutoffBK(1)) * outcomes(1,2);
        AHARV(:,i) = 1*(XbHARV <= cutoffHARV(1)) * outcomes(1,3);
        AMCD(:,i) = 1*(XbMCD <= cutoffMCD(1)) * outcomes(1,4);
        AWEND(:,i) = 1*(XbWEND <= cutoffWEND(1)) * outcomes(1,5);
    elseif i == ng+1
        AAW(:,i) = 1*(XbAW > cutoffAW(ng)) * outcomes(ng+1,1);
        ABK(:,i) = 1*(XbBK > cutoffBK(ng)) * outcomes(ng+1,2);
        AHARV(:,i) = 1*(XbHARV > cutoffHARV(ng)) * outcomes(ng+1,3);
        AMCD(:,i) = 1*(XbMCD > cutoffMCD(ng)) * outcomes(ng+1,4);
        AWEND(:,i) = 1*(XbWEND > cutoffWEND(ng)) * outcomes(ng+1,5);
    else
        AAW(:,i) = 1*(XbAW > cutoffAW(i-1) & XbAW <= cutoffAW(i)) * outcomes(i,1);
        ABK(:,i) = 1*(XbBK > cutoffBK(i-1) & XbBK <= cutoffBK(i)) * outcomes(i,2);
        AHARV(:,i) = 1*(XbHARV > cutoffHARV(i-1) & XbHARV <= cutoffHARV(i)) * outcomes(i,3);
        AMCD(:,i) = 1*(XbMCD > cutoffMCD(i-1) & XbMCD <= cutoffMCD(i)) * outcomes(i,4);
        AWEND(:,i) = 1*(XbWEND > cutoffWEND(i-1) & XbWEND <= cutoffWEND(i)) * outcomes(i,5);
    end
end

% Collapse matrix into single vector.
A = [AAW; ABK; AHARV; AMCD; AWEND];
A = sum(A,2);

end
