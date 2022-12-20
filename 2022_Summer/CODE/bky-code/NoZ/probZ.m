% Dynamic Spillovers in the Retail Industry.        %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.  %
% Returns the probabilities and index values.       %
% November 27, 2012.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fy = probZ(theta, y, X)

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
Xb = X*beta;
N = length(Xb) / 5;
XbAW = Xb(1:N,:);
XbBK = Xb(N+1:2*N,:);
XbHARV = Xb(2*N+1:3*N,:);
XbMCD = Xb(3*N+1:4*N,:);
XbWEND = Xb(4*N+1:5*N,:);

% Get chain specific dependent variables.
yAW = y(1);
yBK = y(2);
yHARV = y(3);
yMCD = y(4);
yWEND = y(5);

% First loop goes through each possible category.
for i = 1:ng+1
    if i == 1
        if (yAW <= grid(1))
            FAW = normcdf(cutoffAW(1) - XbAW);
        end
        if (yBK <= grid(1))
            FBK = normcdf(cutoffBK(1) - XbBK);
        end
        if (yHARV <= grid(1))
            FHARV = normcdf(cutoffHARV(1) - XbHARV);
        end
        if (yMCD <= grid(1))
            FMCD = normcdf(cutoffMCD(1) - XbMCD);
        end
        if (yWEND <= grid(1))
            FWEND = normcdf(cutoffWEND(1) - XbWEND);
        end
    elseif i == ng+1
        if (yAW > grid(ng))
            FAW = 1 - normcdf(cutoffAW(ng) - XbAW);
        end
        if (yBK > grid(ng))
            FBK = 1 - normcdf(cutoffBK(ng) - XbBK);
        end
        if (yHARV > grid(ng))
            FHARV = 1 - normcdf(cutoffHARV(ng) - XbHARV);
        end
        if (yMCD > grid(ng))
            FMCD = 1 - normcdf(cutoffMCD(ng) - XbMCD);
        end
        if (yWEND > grid(ng))
            FWEND = 1 - normcdf(cutoffWEND(ng) - XbWEND);
        end
    else
        if (yAW > grid(i-1) & yAW <= grid(i))
            FAW = normcdf(cutoffAW(i) - XbAW) - normcdf(cutoffAW(i-1) - XbAW);
        end
        if (yBK > grid(i-1) & yBK <= grid(i))
            FBK = normcdf(cutoffBK(i) - XbBK) - normcdf(cutoffBK(i-1) - XbBK);
        end
        if (yHARV > grid(i-1) & yHARV <= grid(i))
            FHARV = normcdf(cutoffHARV(i) - XbHARV) - normcdf(cutoffHARV(i-1) - XbHARV);
        end
        if (yMCD > grid(i-1) & yMCD <= grid(i))
            FMCD = normcdf(cutoffMCD(i) - XbMCD) - normcdf(cutoffMCD(i-1) - XbMCD);
        end
        if (yWEND > grid(i-1) & yWEND <= grid(i))
            FWEND = normcdf(cutoffWEND(i) - XbWEND) - normcdf(cutoffWEND(i-1) - XbWEND);
        end
    end
end

Fy = [FAW; FBK; FHARV; FMCD; FWEND];

end
