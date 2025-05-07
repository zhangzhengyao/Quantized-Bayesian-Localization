function G = jacob(X, x)
% Jacobian matrix computation
% --------------------------------
% G = jacobian(X, x)
% G = Jacobian matrix 
% x = 2D position estimate
% X = matrix for receiver positions
%

[dim,L] = size(X); % L is number of receivers; dim is dimension of space
f_TOA = sqrt(sum((ones(L,1)*x'-X').^2,2));
G = (ones(L,1)*x'- X')./(f_TOA*ones(1,dim));
