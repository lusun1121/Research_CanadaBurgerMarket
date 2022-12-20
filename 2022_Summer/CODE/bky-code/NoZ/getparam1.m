% Dynamic Spillovers in the Retail Industry.            %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.      %
% Parses out first stage parameters.                    %
% July 22, 2013.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [theta, mu, delta, beta, eta] = getparam1(param, specification)

% Constants
numplayers = 5;                         % Number of players
nummarkets = 31;                        % Number of markets
ng = 7;                                 % Number of ordered probit cutoffs

% Ordered probit parameters ($\theta$)
%
% * Coefficients for each firm:
%
%     1. Population
%     2. Income
%     3. Property value
%     4. Grey Cup host
%     5. Smoking regulation
%     6. Minimum wage
%     7. Population * Population
%     8. Population * Income
%     9. Population * Property value
%     10. Population * Minimum wage
%     11. Income * Income
%     12. Income * Property value
%     13. Income * Minimum wage
%     14. Property value * Property value
%     15. Property value * Minimum wage
%     16. Minimum wage * Minimum wage
%     17. Z
%     18. Z * Population
%     19. Z * Income
%     20. Z * Property value
%     21. Z * Minimum wage
%     22. Rival Z
%     23. Rival Z * Population
%     24. Rival Z * Income
%     25. Rival Z * Property value
%     26. Rival Z * Minimum wage
%     27. Lagged N
%
% * Cutoffs for each firm (ng of them).
%
% Overall length of $\theta$: (27+ng) * numplayers.

% Lengths of some internal parameter vectors are fixed
numtheta_total = 27;
nummu = numplayers;
numdelta = numplayers;
numeta = nummarkets;
idx = 0;

% Retrieve parameters depending on specification
if  specification == 1

    % Baseline homogeneous specification (57 parameters).
    % No latent variable, but firm and city fixed effects are included (through Z_imt).

    if (size(param, 1) ~= 50 + ng)
        error('Incorrect number of first stage parameters (specification 3)!')
    end

    % Ordered probit coefficients ($\theta$)
    numtheta_coeff = 16;
    numtheta_zeros = numtheta_total - numtheta_coeff - 1;
    theta_coeff = [ param(1:numtheta_coeff); 1; zeros(numtheta_zeros, 1) ];
    idx = idx + numtheta_coeff;

    % Ordered probit cutoffs ($\theta$)
    numtheta_cutoff = ng;
    theta_cutoff = param(idx+1:idx+numtheta_cutoff);
    idx = idx + numtheta_cutoff;

    % Combine to build theta vector
    theta = [ theta_coeff; theta_cutoff ];
    theta = [ theta; theta; theta; theta; theta; ];

    % Firm-specific drift parameters ($\mu$) (AW = 0)
    mu = [ 0; param(idx+1:idx+numplayers-1) ];
    idx = idx + (numplayers - 1);

    % Firm-specific autoregressive parameters ($\delta$)
    delta = zeros(5,1);

    % Remaining firm-specific coefficients in Z process ($\beta$)
    beta = zeros(6,1);
    beta = [ beta; beta; beta; beta; beta ];

    % City fixed effects in Z process ($\eta)
    eta = [ 0; param(idx+1:idx+nummarkets-1); ];
    idx = idx + (nummarkets - 1);

elseif specification == 2

    % Baseline homogeneous specification with Z interactions (72 parameters).
    % No spillovers, but firm and city fixed effects are included.

    if (size(param, 1) ~= 65 + ng)
        error('Incorrect number of first stage parameters (specification 16)!')
    end

    % Ordered probit coefficients ($\theta$)
    numtheta_coeff = 16;
    numtheta_coeff2 = 9;
    theta_coeff = [ param(1:numtheta_coeff); 1;
                    param(numtheta_coeff+1:numtheta_coeff+numtheta_coeff2); 0 ];
    idx = idx + numtheta_coeff + numtheta_coeff2;

    % Ordered probit cutoffs ($\theta$)
    numtheta_cutoff = ng;
    theta_cutoff = param(idx+1:idx+numtheta_cutoff);
    idx = idx + numtheta_cutoff;

    % Combine to build theta vector
    theta = [ theta_coeff; theta_cutoff ];
    theta = [ theta; theta; theta; theta; theta; ];

    % Firm-specific drift parameters ($\mu$) (AW = 0)
    mu = [ 0; param(idx+1:idx+numplayers-1) ];
    idx = idx + (numplayers - 1);

    % Firm-specific autoregressive parameters ($\delta$)
    delta = abs(param(idx+1:idx+numplayers));
    idx = idx + numplayers;

    % Remaining firm-specific coefficients in Z process ($\beta$)
    numbeta = 1;
    tmp = param(idx+1);
    beta = [ 0; 0; 0; 0; 0; abs(tmp) ];
    beta = [ beta; beta; beta; beta; beta ];
    idx = idx + numbeta;

    % City fixed effects in Z process ($\eta)
    eta = [ 0; param(idx+1:idx+nummarkets-1); ];
    idx = idx + (nummarkets - 1);

elseif specification == 3

    % Heterogeneous cubic N specification (91 parameters).
    % Includes spillovers, rival spillovers, spillover interactions, and
    % firm and city fixed effects.

    if (size(param, 1) ~= 84 + ng)
        error('Incorrect number of first stage parameters (specification 14)!')
    end

    % Ordered probit coefficients ($\theta$)
    numtheta_coeff = 16;
    numtheta_coeff2 = 9;
    theta_coeff = [ param(1:numtheta_coeff); 1;
                    param(numtheta_coeff+1:numtheta_coeff+numtheta_coeff2); 0 ];
    idx = idx + numtheta_coeff + numtheta_coeff2;

    % Ordered probit cutoffs ($\theta$)
    numtheta_cutoff = ng;
    theta_cutoff = param(idx+1:idx+numtheta_cutoff);
    idx = idx + numtheta_cutoff;

    % Combine to build theta vector
    theta = [ theta_coeff; theta_cutoff ];
    theta = [ theta; theta; theta; theta; theta; ];

    % Firm-specific drift parameters ($\mu$) (AW = 0)
    mu = [ 0; param(idx+1:idx+numplayers-1) ];
    idx = idx + (numplayers - 1);

    % Firm-specific autoregressive parameters ($\delta$)
    delta = abs(param(idx+1:idx+numplayers));
    idx = idx + numplayers;

    % Remaining firm-specific coefficients in Z process ($\beta$)
    numbeta = 4 * numplayers;
    tmp = param(idx+1:idx+numbeta);
    beta = [ 0; 0; tmp(1);  tmp(2);  tmp(3);  abs(tmp(4));
             0; 0; tmp(5);  tmp(6);  tmp(7);  abs(tmp(8));
             0; 0; tmp(9);  tmp(10); tmp(11); abs(tmp(12));
             0; 0; tmp(13); tmp(14); tmp(15); abs(tmp(16));
             0; 0; tmp(17); tmp(18); tmp(19); abs(tmp(20)); ];
    idx = idx + numbeta;

    % City fixed effects in Z process ($\eta)
    eta = [ 0; param(idx+1:idx+nummarkets-1); ];
    idx = idx + (nummarkets - 1);

end
