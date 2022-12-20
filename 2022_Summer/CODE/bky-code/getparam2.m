% Dynamic Spillovers in the Retail Industry.            %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.      %
% Parses out second stage parameters.                   %
% September 23, 2014.                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function thetaS = getparam2(param, specification)

% Constants
numplayers = 5;                         % Number of players

% Profit function parameters:
%
% * $\theta_1$: X's in levels (6 variables)
% * $\theta_2$: X's interacted with N (6 variables)
% * $\gamma$: coefficient on spillovers
% * $\theta_3$: competitive effect (7 variables: level, log, AW, BK, HARV, MCD, WEND)
% * $\theta_4$: own stores (3 variables: level, square, log)
% * $\psi_1$: entry cost
% * $\psi_2$: building cost
% * $\psi_3$: scrap value
% * $1$: coefficient on error
%
% Total number of parameters: 22 * numplayers = 110

% Retrieve parameters depending on specification
if specification == 1

    % Homogeneous specification (11 parameters).
    % Parameters are identical across firms.
    % Competitive and own store effects are in levels (no square).
    % Separate entry cost.  No latent variable.

    if size(param, 1) ~= 11
        error(['Incorrect number of second stage parameters!'])
    end

    theta1 = param(1:6);
    theta2 = zeros(6,1);
    gamma = 1; % Coefficient on Z, which is zero or a sum of fixed effects
    theta3 = [ param(7); 0; 0; 0; 0; 0; 0 ];
    theta4 = [ param(8); 0; 0 ];
    psi1 = param(9);
    psi2 = param(10);
    psi3 = param(11);
    tmp = [ theta1; theta2; gamma; theta3; theta4; psi1; psi2; psi3; 1 ];
    thetaS = [ tmp; tmp; tmp; tmp; tmp; ];

elseif specification == 2 || specification == 3

    % Homogeneous specification (12 parameters).
    % Parameters are identical across firms.
    % Competitive and own store effects are in levels (no square).
    % Separate entry cost.

    if size(param, 1) ~= 12
        error(['Incorrect number of second stage parameters!'])
    end

    theta1 = param(1:6);
    theta2 = zeros(6,1);
    gamma = param(7);
    theta3 = [ param(8); 0; 0; 0; 0; 0; 0 ];
    theta4 = [ param(9); 0; 0 ];
    psi1 = param(10);
    psi2 = param(11);
    psi3 = param(12);
    tmp = [ theta1; theta2; gamma; theta3; theta4; psi1; psi2; psi3; 1 ];
    thetaS = [ tmp; tmp; tmp; tmp; tmp; ];

end
