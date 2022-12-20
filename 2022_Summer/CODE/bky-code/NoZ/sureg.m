% Dynamic Spillovers in the Retail Industry.        %
% By Jason Blevins, Ahmed Khwaja, and Nathan Yang.  %
% Seemingly unrelated regression estimation.        %
% September 20, 2012.                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results, sigma] = sureg(X,lagX)

% Generate column of ones for constant.
nobs = size(X,1);
neqs = size(X,2);

% Vector of ones for the constant.
o = ones(nobs,1);

% Reshape the dependent variables, and its lag.
y = reshape(X,nobs*neqs,1);
M = [lagX(:,1:neqs),o];
zM = zeros(size(M,1),size(M,2));
A = [M, zM, zM, zM; zM, M, zM, zM; zM, zM, M, zM; zM, zM, zM, M];

% SUR regression
bOLS = (A'*A)^(-1)*A'*y;

% Calculate the residuals.
u = y - A * bOLS;

% Compute SUR residuals.
emat = zeros(nobs,neqs);

% Parse out the residuals based on equation.
for i = 1:neqs
    emat(:,i) = u(nobs*(i-1)+1:nobs*i);
end

% Compute SUR sigma matrix.
sigma = zeros(neqs,neqs);

for i=1:neqs;
 for j=i:neqs;
        sigma(i,j) = (emat(:,i)'*emat(:,j))/nobs;
    if j > 1;
        sigma(j,i) = sigma(i,j);
    end;
 end;
end;

sigmai = inv(sigma);

nx = neqs * size(A,2);

% compute sur var-cov matrix
xx = zeros(nx,nx);

for i=1:neqs;
 for j=1:neqs;
  xx = kron(sigmai,(M'*M));
 end;
end;

% Find inverse of var-cov matrix for inference.
xxi = inv(xx);

% Compute standard errors from the SUR regression.
vcov = diag(xxi);
std = vcov.^(1/2);

% Compile the results for the parameter estimates and standard errors.
results = [bOLS, std];

% Calculate the variance matrix for the errors.
sigma = cov(emat);

end
