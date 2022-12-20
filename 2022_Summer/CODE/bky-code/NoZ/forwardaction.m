% Dynamic Spillovers in the Retail Industry.        %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.  %
% Simulates the choices for each simulation.        %
% November 14, 2012.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, Xb] = forwardaction(theta,X,u)

% Construct a grid of cutoffs.
grid = [-5 -2 -1 0 1 2 5];
ng = length(grid);

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
for i = 1:ng
    if i == 1
        AAW(:,i) = 1*(XbAW <= cutoffAW(1)) * grid(1);
        ABK(:,i) = 1*(XbBK <= cutoffBK(1)) * grid(1);
        AHARV(:,i) = 1*(XbHARV <= cutoffHARV(1)) * grid(1);
        AMCD(:,i) = 1*(XbMCD <= cutoffMCD(1)) * grid(1);
        AWEND(:,i) = 1*(XbWEND <= cutoffWEND(1)) * grid(1);
    elseif i == ng
        AAW(:,i) = 1*(XbAW > cutoffAW(ng)) * grid(ng);
        ABK(:,i) = 1*(XbBK > cutoffBK(ng)) * grid(ng);
        AHARV(:,i) = 1*(XbHARV > cutoffHARV(ng)) * grid(ng);
        AMCD(:,i) = 1*(XbMCD > cutoffMCD(ng)) * grid(ng);
        AWEND(:,i) = 1*(XbWEND > cutoffWEND(ng)) * grid(ng);
    else
        AAW(:,i) = 1*(XbAW > cutoffAW(i-1) & XbAW <= cutoffAW(i)) * grid(i);
        ABK(:,i) = 1*(XbBK > cutoffBK(i-1) & XbBK <= cutoffBK(i)) * grid(i);
        AHARV(:,i) = 1*(XbHARV > cutoffHARV(i-1) & XbHARV <= cutoffHARV(i)) * grid(i);
        AMCD(:,i) = 1*(XbMCD > cutoffMCD(i-1) & XbMCD <= cutoffMCD(i)) * grid(i);
        AWEND(:,i) = 1*(XbWEND > cutoffWEND(i-1) & XbWEND <= cutoffWEND(i)) * grid(i);
    end
end

% Collapse matrix into single vector.
A = [AAW; ABK; AHARV; AMCD; AWEND];
A = sum(A,2);

end
