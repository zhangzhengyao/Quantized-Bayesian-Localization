function H = hessian_ml(X,x,r,sigma2)
% ML Hessian matrix computation
% --------------------------------
% H = hessian_ml(X,x,r,sigma2)
% H = Hessian matrix 
% X = matrix for receiver positions
% x = 2D position estimate
% r = TOA measurement vector
% sigma2 = noise variance vector
%
L = size(X,2); % number of receivers

t1 = 0;
t2 = 0;
t3 =0;
ds = sum((x*ones(1,L)-X).^2,1);
ds = ds';
for i=1:L
    t1 = t1 + (1/sigma2(i))*((x(1)-X(1,i))^2/ds(i)-(r(i)-ds(i)^(0.5))*(x(2)-X(2,i))^2/ds(i)^(1.5));
    t2 = t2 + (1/sigma2(i))*((x(2)-X(2,i))^2/ds(i)-(r(i)-ds(i)^(0.5))*(x(1)-X(1,i))^2/ds(i)^(1.5));
    t3 = t3 + (1/sigma2(i))*(r(i)*(x(1)-X(1,i))*(x(2)-X(2,i))/ds(i)^(1.5));
end
H=2.*[t1 t3;
      t3 t2];