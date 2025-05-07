function x = gn_ml(X,r,iter,sigma2)
% Maximum Likelihood
% Local search scheme: Gauss-Netwon algorithm
% --------------------------------
% x = gn_ml(X,r,iter,sigma2);
% x = 2D position estimate
% X = % matrix for receiver positions
% r = TOA measurement vector
% iter = number of iterations
% sigma2 = noise variance vector
% 
L = size(X,2); % number of receivers
x = wlls(X,r,sigma2);
for i = 1:iter
    G = jacob(X, x);
    f_TOA = sqrt(sum((ones(L,1)*x'-X').^2,2));
    C_inv = diag(1./sigma2);
    x = x+inv(G'*C_inv*G)*G'*C_inv*(r-f_TOA);
end