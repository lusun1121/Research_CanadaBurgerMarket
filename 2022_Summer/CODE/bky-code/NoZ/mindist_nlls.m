% Dynamic Spillovers in the Retail Industry.                    %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.              %
% Calculates the BBL minimum distance objective function.       %
% August 7, 2013.                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F, J] = mindist_nlls(param, APIstar, APIperturb, secondspec)

numplayers = 5;
nummarkets = 31;
numineq = size(APIstar,1);
numparam = size(param,1);

thetaS = getparam2(param, secondspec);

% Difference in coefficients betweeen equilibrium and perturbed policies
APIdiff = APIstar - APIperturb;
Vdiff = APIdiff*thetaS;

% Initialize criterion, gradient, and Hessian.
F = zeros(numineq, 1);

% Calculate nonlinear least square component functions
Ibind = Vdiff < 0;
F = Ibind .* Vdiff;

% Calculate Jacobian (if requested)
if nargout > 1
    J = zeros(numineq, numparam);
    for j = 1:numineq
        if Ibind(j) == 1
            J(j,:) = APIdiff(j,:)*PS;
        end
    end
end

end
