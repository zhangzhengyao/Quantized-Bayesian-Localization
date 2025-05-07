function g = grad_ml(X,x,r,sigma2)
% ML gradient computation
% --------------------------------
% g = grad_ml(X,x,r);
% g = gradient vector 
% X = matrix for receiver positions
% x = 2D position estimate
% r = TOA measurement vector
% sigma2 = noise variance vector
%
L = size(X,2); % number of receivers
t1 = 0;
t2 = 0;
ds = sum((x*ones(1,L)-X).^2,1);
ds = ds';
for i=1:L
    t1 = t1 + (1/sigma2(i))*(r(i)-ds(i)^(0.5))*(x(1)-X(1,i))/ds(i)^(0.5);
    t2 = t2 + (1/sigma2(i))*(r(i)-ds(i)^(0.5))*(x(2)-X(2,i))/ds(i)^(0.5);
end
g=-2.*[t1; t2];